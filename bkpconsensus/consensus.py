import fire
import pandas as pd
from bkpconsensus.breakpoint_db import BreakpointDatabase
from single_cell.utils import csvutils
from bkpconsensus.vcf_sv_parser import SvVcfData


def read_destruct(destruct_calls):
    df = csvutils.CsvInput(destruct_calls).read_csv()
    df['breakpoint_id'] = df['prediction_id']
    return df


def read_lumpycsv(lumpy_filename):
    df = csvutils.CsvInput(lumpy_filename).read_csv()
    df = df.rename(columns={
        'chrom1': 'chromosome_1',
        'strand1': 'strand_1',
        'chrom2': 'chromosome_2',
        'strand2': 'strand_2',
    })
    df['position_1'] = 0.5 * (df['confidence_interval_start1'] + df['confidence_interval_end1'])
    df['position_2'] = 0.5 * (df['confidence_interval_start2'] + df['confidence_interval_end2'])
    return df


def read_consensus(filename):
    try:
        data = pd.read_csv(filename, sep=',', compression=None, converters={'chromosome_1': str, 'chromosome_2': str})
        data['breakpoint_id'] = data['prediction_id']
    except UnicodeDecodeError:
        data = pd.read_csv(filename, sep=',', compression='gzip', converters={'chromosome_1': str, 'chromosome_2': str})
    return data


def check_common(x, df_db, calls):
    val = df_db.query(x, extend=500)

    val = sorted(val)

    if len(val) == 1:
        return

    if val[0] not in calls:
        calls[val[0]] = set()

    for v in val[1:]:
        calls[val[0]].add(v)


def get_common_calls(df, df_db):
    calls = {}

    for i, row in df.iterrows():
        check_common(row, df_db, calls)

    new_groups = {}
    for i, (key, vals) in enumerate(calls.items()):
        new_groups[key] = i
        for val in vals:
            new_groups[val] = i

    return new_groups


def consensusmulti(destruct_calls, lumpy_calls, svaba_calls, gridss_calls, consensus_calls):
    allcalls = [
        read_destruct(destruct_calls),
        SvVcfData(lumpy_calls).as_data_frame(),
        SvVcfData(svaba_calls).as_data_frame(),
        SvVcfData(gridss_calls).as_data_frame()
    ]

    allcalls = pd.concat(allcalls)

    allcalls_db = BreakpointDatabase(allcalls)

    groups = get_common_calls(allcalls, allcalls_db)

    allcalls['grouped_breakpoint_id'] = allcalls['breakpoint_id'].apply(lambda x: groups.get(x, float("nan")))

    allcalls = allcalls[~ pd.isnull(allcalls.grouped_breakpoint_id)]

    allcalls = allcalls.groupby('grouped_breakpoint_id')

    outdata = []
    for _, brkgrp in allcalls:

        # filter multiple calls by same tool in the window
        # without confirmation from another tool
        if len(brkgrp.caller.unique()) == 1:
            continue

        brkgrp['caller'] = ','.join(list(brkgrp['caller']))
        brkgrp = brkgrp[:1]

        outdata.append(brkgrp)

    outdata = pd.concat(outdata)

    outdata.to_csv(consensus_calls, index=False)


def read_sv(filename, type_):
    if type_ == 'destruct':
        data = read_destruct(filename)
    elif type_ == 'lumpycsv':
        data = read_lumpycsv(filename)
    elif type_ == 'consensus':
        data = read_consensus(filename)
    elif type_ in ('lumpy', 'svaba', 'gridss'):
        data = SvVcfData(filename).as_data_frame()
    else:
        raise ValueError(f'unrecognized type {type_}')
    return data


def consensus(filename1, type1, filename2, type2, out_filename, min_dist=200):
    data1 = read_sv(filename1, type1)
    data2 = read_sv(filename2, type2)

    if 'breakpoint_id' not in data1:
        raise ValueError('breakpoint_id not in data1')

    if 'breakpoint_id' not in data2:
        raise ValueError('breakpoint_id not in data2')

    db1 = BreakpointDatabase(data1)

    results = []
    for idx, row in data2.iterrows():
        for match_id in db1.query(row, min_dist):
            results.append({
                'breakpoint_id_1': match_id,
                'breakpoint_id_2': row['breakpoint_id']})
    results = pd.DataFrame(results)

    if results.empty:
        results = pd.DataFrame(columns=['breakpoint_id_1', 'breakpoint_id_2'])

    results.to_csv(out_filename, index=False)


def main():
    fire.Fire(consensus)




import fire
import pandas as pd
from bkpconsensus.breakpoint_db import BreakpointDatabase
from single_cell.utils import csvutils
from bkpconsensus.vcf_sv_parser import SvVcfData


def read_destruct(destruct_calls):
    df = csvutils.CsvInput(destruct_calls).read_csv()
    destruct_cols = ['prediction_id', 'chromosome_1', 'position_1', 'strand_1', 'chromosome_2', 'position_2', 'strand_2', 'type']
    df = df[destruct_cols]
    df['breakpoint_id'] = df['prediction_id']
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


def consensus(filename1, type1, filename2, type2, out_filename, min_dist=200):
    if type1 == 'destruct':
        data1 = read_destruct(filename1)
    elif type1 == 'consensus':
        data1 = read_consensus(filename1)
    elif type1 in ('lumpy', 'svaba', 'gridss'):
        data1 = SvVcfData(filename1).as_data_frame()
    else:
        raise ValueError(f'unrecognized type {type1}')

    if 'breakpoint_id' not in data1:
        raise ValueError('breakpoint_id not in data1')

    if type2 == 'destruct':
        data2 = read_destruct(filename2)
    elif type2 == 'consensus':
        data2 = read_consensus(filename2)
    elif type2 in ('lumpy', 'svaba', 'gridss'):
        data2 = SvVcfData(filename2).as_data_frame()
    else:
        raise ValueError(f'unrecognized type {type2}')

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




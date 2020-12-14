import pandas as pd
from breakpoint_db import BreakpointDatabase
from single_cell.utils import csvutils
from vcf_sv_parser import SvVcfData


def read_destruct(destruct_calls):
    df = csvutils.CsvInput(destruct_calls).read_csv()
    df = df[
        ['prediction_id', 'chromosome_1', 'position_1', 'strand_1', 'chromosome_2', 'position_2', 'strand_2', 'type']]
    df['breakpoint_id'] = df.prediction_id
    del df['prediction_id']
    df['caller'] = 'destruct'

    df['breakpoint_id'] = df['breakpoint_id'].astype(str) + '_' + df['caller']
    return df


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


def consensus(destruct_calls, lumpy_calls, svaba_calls, gridss_calls, consensus_calls):
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

consensus(
    "../wgs/output/breakpoints/SA1256PP/SA1256PP_destruct_breakpoints.csv.gz",
    '../wgs/output/breakpoints/SA1256PP/SA1256PP_lumpy.vcf',
    '../svaba/output/out.svaba.somatic.sv.vcf.gz',
    '../gridss/calls.vcf.gz',
    'consensus.csv.gz'
)

def main():
    print('main')




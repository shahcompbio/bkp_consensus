import scgenome
import scgenome.db.qc_from_files
import scgenome.cnfilter
import pandas as pd
import scipy.stats
import scgenome.loaders.breakpoint
import statistics
import scgenome.breakpointdata
import numpy as np 
import single_cell.utils.csvutils
from single_cell.utils.csvutils import IrregularCsvInput, CsvOutput
from vcf_sv_parser import SvVcfData


def load_qc_data(
        sample_id, annotation_metrics, hmmcopy_reads, hmmcopy_segs,
        hmmcopy_metrics, alignment_metrics, gc_metrics
):
    results_tables_new = scgenome.db.qc_from_files.get_qc_data_from_filenames(
        [annotation_metrics],
        [hmmcopy_reads],
        [hmmcopy_segs],
        [hmmcopy_metrics], [alignment_metrics], [gc_metrics],
        sample_ids=[sample_id], additional_hmmcopy_reads_cols=None
    )
    cn_data = results_tables_new['hmmcopy_reads']
    metrics_data = results_tables_new['annotation_metrics']
    metrics_data = scgenome.cnfilter.calculate_filter_metrics(
        metrics_data,
        cn_data,
    )
    filtered_cells = metrics_data.loc[
        metrics_data['filter_is_s_phase'] &
        metrics_data['filter_quality'],
        ['cell_id']]
    cn_data_filt = cn_data.merge(filtered_cells[['cell_id']].drop_duplicates())
    return metrics_data, cn_data_filt

def get_cna_changepoints(cna_data):
    cn_data_summary = cna_data.groupby(['chr', 'start', 'end'])
    cn_data_summary = cn_data_summary[['copy', 'state']]
    cn_data_summary = cn_data_summary.aggregate(
        {'copy': statistics.median, 'state': statistics.median}
    )
    cn_data_summary = cn_data_summary.sort_values(by=["chr", "start"])
    cn_data_summary["cn_change"]=cn_data_summary.groupby('chr').state.transform(pd.Series.diff)
    cn_data_summary.reset_index(inplace=True)
    #get the places where theres a cn change and record both hose locations 
    #and the locations before them to include both indexes involved in the cn change
    change_data = cn_data_summary[(cn_data_summary.cn_change!=0) & (cn_data_summary.cn_change!=np.NaN)].cn_change
    change_data.index = change_data.index-1
    cn_data_summary.update(change_data)
    cn_data_summary = cn_data_summary[cn_data_summary.cn_change != 0 ]
    cn_data_summary.reset_index(inplace=True)
    return cn_data_summary

def anno_bkp(row, cn, side="1"):
    if side == "1":
        chrom = str(row.chromosome_1)
        pos = row.position_1
    else:
        chrom = str(row.chromosome_2)
        pos = row.position_2        
    pos = int(pos)
    matches = cn[(cn.start < pos) 
        & (cn.end > pos) 
        & (cn.chr == chrom)]
    matches = matches.reset_index()
    if matches.empty:
        row["chr_pos_{}".format(side)] = np.NaN
        row["start_pos_{}".format(side)] = np.NaN
        row["end_pos_{}".format(side)] = np.NaN
        row["state_pos_{}".format(side)] = np.NaN
        row["cn_change_pos_{}".format(side)] = np.NaN
    else:
        row["chr_pos_{}".format(side)] = matches.chr.tolist()[0]
        row["start_pos_{}".format(side)] = matches.start.tolist()[0]
        row["end_pos_{}".format(side)] = matches.end.tolist()[0]
        row["state_pos_{}".format(side)] = matches.state.tolist()[0]
        row["cn_change_pos_{}".format(side)] = matches.cn_change.tolist()[0]
    row["pos_{}_in_cn_change_region".format(side)] = not matches.empty
    return row

def match_changepoints(changepoints, sv):
    sv = sv.apply(lambda row: anno_bkp(row, changepoints, side="2"), axis=1)
    sv = sv.reset_index()
    sv = sv.apply(lambda row: anno_bkp(row, changepoints), axis=1)
    return sv

def load_breakpoint_data(
        sample_id, library_id, breakpoint_annotation, breakpoint_count,
        lumpy=False, filter_data=False
):
    breakpoint_results = scgenome.loaders.breakpoint.load_breakpoint_data_from_files(
        [breakpoint_annotation],
        [breakpoint_count],
        lumpy=lumpy
    )
    breakpoint_data = breakpoint_results['breakpoint_data']
    breakpoint_data["library_id"] = library_id
    breakpoint_data["sample_id"] = sample_id
    breakpoint_count_data = breakpoint_results['breakpoint_count_data']
    breakpoint_count_data["library_id"] = library_id
    breakpoint_count_data["sample_id"] = sample_id
    breakpoint_data, breakpoint_count_data = scgenome.breakpointdata.annotate_breakpoint_data(
        breakpoint_data,
        breakpoint_count_data,
        is_lumpy=lumpy
    )
    if filter_data:
        breakpoint_data, breakpoint_count_data = scgenome.breakpointdata.filter_breakpoint_data(
            breakpoint_data, breakpoint_count_data,
        )
    return breakpoint_data

def get_concordance(sv, changepoints):
    print(sv)
    sv = single_cell.utils.csvutils.CsvInput(sv).read_csv()
    print("Total: {}".format(len(sv)))
    sv = match_changepoints(changepoints,
            sv[["chromosome_1", "chromosome_2", "position_1", "position_2"]]
        )
    n_concordant = sv[(sv.pos_1_in_cn_change_region == True) & (sv.pos_2_in_cn_change_region  == True)]
    print("HMM-concordant: {}".format(len(n_concordant)))

def write_metadata(infile, dtypes):
    csvinput = IrregularCsvInput(infile, dtypes)
    csvoutput = CsvOutput(
        infile, csvinput.dtypes, header=csvinput.header,
        columns=csvinput.columns
    )
    csvoutput.write_yaml()


sample_id = "SA1256PP"
library_id = "A108838B"
annotation_metrics= "/juno/work/shah/tantalus/SC-3803/results/annotation/A108838B_metrics.csv.gz"
hmmcopy_reads = "/juno/work/shah/tantalus/SC-3803/results/hmmcopy/A108838B_reads.csv.gz"
hmmcopy_segs = "/juno/work/shah/tantalus/SC-3803/results/hmmcopy/A108838B_segments.csv.gz"
hmmcopy_metrics = "/juno/work/shah/tantalus/SC-3803/results/hmmcopy/A108838B_hmmcopy_metrics.csv.gz"
alignment_metrics = "/juno/work/shah/tantalus/SC-3803/results/align/A108838B_alignment_metrics.csv.gz"
gc_metrics = "/juno/work/shah/tantalus/SC-3803/results/align/A108838B_gc_metrics.csv.gz"

dtypes={
    'chromosome_1': int,
    'position_1': int,
    'chromosome_2': int,
    'position_2': int,
    'strand_1': int,
    'strand_2': int,
    'type': str,
    'caller': str,
    'breakpoint_id': str,
    'grouped_breakpoint_id': str
}

metrics_data, cn_data_filt = load_qc_data(
        sample_id, annotation_metrics, hmmcopy_reads,
        hmmcopy_segs, hmmcopy_metrics, alignment_metrics, gc_metrics)

changepoints = get_cna_changepoints(cn_data_filt)

get_concordance("../wgs/output/breakpoints/SA1256PP/SA1256PP_destruct_breakpoints.csv.gz", changepoints)


lumpy_data = "../wgs/output/breakpoints/SA1256PP/SA1256PP_lumpy.vcf"
df = SvVcfData(lumpy_data).as_data_frame()
df.to_csv("lumpy_parsed.csv.gz", index=False)
write_metadata("lumpy_parsed.csv.gz", dtypes)
get_concordance("lumpy_parsed.csv.gz",changepoints)


svaba_data = "../svaba/output/out.svaba.somatic.sv.vcf.gz"
df = SvVcfData(svaba_data).as_data_frame()
df.to_csv("svaba_parsed.csv.gz", index=False)
write_metadata("svaba_parsed.csv.gz", dtypes)
get_concordance("svaba_parsed.csv.gz",changepoints)


gridss_data = "../gridss/calls.vcf.gz"
df = SvVcfData(gridss_data).as_data_frame()
df.to_csv("gridss_parsed.csv.gz", index=False)
write_metadata("gridss_parsed.csv.gz", dtypes)
get_concordance("gridss_parsed.csv.gz",changepoints)


consensus_calls = "../consensus.csv.gz"
write_metadata("consensus.csv.gz", dtypes)
get_concordance("consensus.csv.gz", changepoints)

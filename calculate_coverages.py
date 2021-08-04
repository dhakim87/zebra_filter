import pandas as pd
from collections import defaultdict
import re
from sys import argv
from os import path
import click
from glob import glob
import gzip
import lzma
from cover import SortedRangeList
from zebra_pdf import compute_probability_coverage_less_than

@click.command()
@click.option('-i',"--input", required=True, help="Input: Directory of sam files (files must end in .sam).")
@click.option('-o',"--output", required=True, help='Output: file name for list of coverages.')
@click.option('-d',"--database", default="databases/WoL/metadata.tsv", help='WoL genome metadata file.',show_default=True)

def calculate_coverages(input, output, database):
    ###################################
    #Calculate coverage of each contig#
    ###################################
    gotu_dict = defaultdict(SortedRangeList)
    gotu_total_reads_dict = defaultdict(int)
    gotu_total_length_dict = defaultdict(int)
    gotu_total_pcr_dups_dict = defaultdict(int)
    file_list = glob(input + "/*.sam")
    file_list_gz = glob(input + "/*.sam.gz")
    file_list_xz = glob(input + "/*.sam.xz")

    file_list = file_list + file_list_gz + file_list_xz
    for samfile in file_list:
        open_sam_file = None
        if samfile.endswith(".sam"):
            open_sam_file = open(samfile.strip(), 'r')
        elif samfile.endswith(".sam.gz"):
            open_sam_file = gzip.open(
                samfile.strip(),
                mode='rt',
                encoding='utf-8')
        elif samfile.endswith(".sam.xz"):
            open_sam_file = lzma.open(
                samfile.strip(),
                mode='rt',
                encoding='utf-8')
        else:
            raise IOError("Unrecognized file extension on '%s'." % samfile)

        pcr_dups_dict = defaultdict(dict)
        with open_sam_file:
            for line in open_sam_file:
                # Get values for contig, location, and length
                linesplit= line.split()
                gotu = linesplit[2]
                location = int(linesplit[3])
                # Get sum of lengths in CIGAR string. Counting deletions as alignment because they should be small
                length_string = linesplit[5]
                length = sum([int(x) for x in re.split("[a-zA-Z]",length_string) if x])
                # Add range to contig_dict
                gotu_dict[gotu].add_range(location, location + length - 1)
                gotu_total_reads_dict[gotu] += 1
                gotu_total_length_dict[gotu] += length
                # Count exact duplicates
                gotu_dups = pcr_dups_dict[gotu]
                if (location,length) not in gotu_dups:
                    gotu_dups[(location, length)] = 1
                else:
                    gotu_dups[(location, length)] += 1

        for gotu_key in pcr_dups_dict:
            all_reads = sum(pcr_dups_dict[gotu_key].values())
            unique_reads = len(pcr_dups_dict[gotu_key])
            gotu_total_pcr_dups_dict[gotu_key] += all_reads - unique_reads


    ###################################
    # Get information from database                #
    ###################################
    md = pd.read_table(database).loc[:,["#genome","total_length","unique_name"]]
    md.columns = ["gotu","total_length","strain"]
    md = md.set_index("gotu")

    #####################
    # Calculate coverages#
    #####################

    num_reads = []
    mean_read_length = []
    estimated_pcr_dups = []
    for k in gotu_dict.keys():
        num_read = gotu_total_reads_dict[k]
        estimated_pcr_dup = gotu_total_pcr_dups_dict[k]
        total_read_length = gotu_total_length_dict[k]
        mean_read_len = 0
        if total_read_length > 0:
            mean_read_len = total_read_length / num_read
        num_reads.append(num_read)
        estimated_pcr_dups.append(estimated_pcr_dup)
        mean_read_length.append(mean_read_len)

    # Make dataframe from dictionary of coverages of each contig
    cov = pd.DataFrame(
        {
            "gotu": list(gotu_dict.keys()),
            "covered_length": [x.compute_length() for x in gotu_dict.values()],
            "num_reads": num_reads,
            "estimated_pcr_dups": estimated_pcr_dups,
            "mean_read_len": mean_read_length
        }
    )
    cov = cov.set_index("gotu")
    cov = cov.sort_values("covered_length", ascending=False)
    # Add genome metadata
    cov = cov.join(md, how="left")
    # Calculate coverage percent
    cov["coverage_ratio"] = cov.apply(func= lambda x : x["covered_length"]/x["total_length"], axis=1)
    cov["p_coverage"] = cov.apply(func = lambda x : compute_probability_coverage_less_than(
        x["coverage_ratio"], x["total_length"], x["mean_read_len"], x["num_reads"] - x["estimated_pcr_dumps"]), axis=1)
    cov = cov.loc[:,["covered_length","total_length","coverage_ratio","strain", "num_reads", "mean_read_len", "p_coverage"]]

    ##############
    # Write output #
    ##############
    cov.to_csv(output, sep='\t')


if __name__ == "__main__":
    calculate_coverages()

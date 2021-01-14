#!/usr/bin/env python3

import os
import sys
import logging

from collections import deque

from .maths import pvalue_adjust

import numpy
import pandas

from scipy.stats import binom_test, chi2_contingency

def _closest(lst, K):
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]


def _f(df):
    return df[df['total_coverage'] == numpy.max(df.total_coverage)]


def _process_ase(genes_df, file_list, min_reads_gene=10):
    # mypath = (data_dir)
    # files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    # for file in files:
    for f in file_list: # type: str
        # if (file.startswith('.')):
        #     continue
        # print("Processing ", f)
        logging.info("Processing %s", f)
        hap_counts_tr = pandas.read_csv(f, sep='\t')  # files should be tab-delimited
        hap_counts_tr['transcript'] = hap_counts_tr['transcript'].astype(str)
        hap_counts_ge = hap_counts_tr.merge(genes_df, on=['transcript', 'transcript'])
        hap_counts_ge['combine'] = hap_counts_ge['gene'].map(str) + \
            "_" + \
            hap_counts_ge['contig'].map(str) + \
            "_" + \
            hap_counts_ge['position'].map(str) + \
            "_" + \
            hap_counts_ge['refAllele'].map(str) + \
            "_" + \
            hap_counts_ge['altAllele'].map(str)
        gene_var_ids = hap_counts_ge['combine'].unique()
        data = deque()
        for gene_var in gene_var_ids:
            hap_count_subset = hap_counts_ge[(hap_counts_ge['combine'] == gene_var)]
            gene_id = gene_var.split('_')[0]
            var_id = gene_var.split('_', 1)[1]
            total_coverage = hap_count_subset['refCount'].sum() + hap_count_subset['altCount'].sum()
            if int(total_coverage) < int(min_reads_gene):
                logging.warning("Not enough reads for %s and %s", gene_id, var_id)
                continue
            real_ref_counts = hap_count_subset['refCount'].sum()
            real_alt_counts = hap_count_subset['altCount'].sum()
            # binom_p = stats.binom_test(real_ref_counts, n=total_coverage, p=0.5)
            binom_p = binom_test(real_ref_counts, n=total_coverage, p=0.5)
            logafc = round(numpy.log2(real_alt_counts / real_ref_counts), 8)
            data.append([gene_id, var_id, total_coverage, logafc, binom_p])
            output_power = pandas.DataFrame({
                'gene_id': [item[0] for item in data], 'variant_id': [item[1] for item in data],
                'total_coverage': [item[2] for item in data], 'logafc': [item[3] for item in data],
                'p_value': [item[4] for item in data]
            })
            # select the SNP with the highest sum count in case of multiple SNP for a gene
            output_power_filtered = (output_power.groupby("gene_id").apply(_f))
            output_power_filtered_indexed = output_power_filtered.set_index('gene_id')
            # output_power_filtered_indexed['qvalue'] = fdr_adjust(output_power_filtered_indexed['p_value'].tolist())
            output_power_filtered_indexed['qvalue'] = pvalue_adjust(pvalues=output_power_filtered_indexed['p_value'].tolist())
            # output_power_filtered_indexed.to_csv(output_dir + "processed_ase_" + file, sep='\t')
            return output_power_filtered_indexed


def _process_asts(genes_df, file_list, min_reads=10):
    for f in file_list:
        logging.info("Processing %s", f)
        hap_counts_tr = pandas.read_csv(f, sep='\t')  # files should be tab-delimited
        hap_counts_tr['transcript'] = hap_counts_tr['transcript'].astype(str)
        hap_counts_ge = hap_counts_tr.merge(genes_df, on=['transcript', 'transcript'])
        hap_counts = hap_counts_ge[(hap_counts_ge['altCount'] >= int(min_reads)) | (hap_counts_ge['refCount'] >= int(min_reads))]
        hap_counts['combine'] = hap_counts['gene'].map(str) + \
            "_" + \
            hap_counts['contig'].map(str) + \
            "_" + \
            hap_counts['position'].map(str) + \
            "_" + \
            hap_counts['refAllele'].map(str) + \
            "_" + \
            hap_counts['altAllele'].map(str)
        gene_var_ids = hap_counts['combine'].unique()
        # print(gene_var_ids)
        data = deque()
        for gene_var in gene_var_ids:
            hap_count_subset = hap_counts[(hap_counts['combine'] == gene_var)]
            number_of_isoforms = len(hap_count_subset)
            gene_id = gene_var.split('_')[0]
            var_id = gene_var.split('_', 1)[1]
            if (number_of_isoforms == 1):
                # print("only one isoform found for", gene_id, " and ", var_id)
                logging.warning("Only one isoform found for %s and %s", gene_id, var_id)
                continue
            if ((hap_count_subset['altCount'] == 0).all() or (hap_count_subset['refCount'] == 0).all()):
                # print("alt or ref have 0 counts for", gene_id, " and ", var_id)
                logging.warning("alt or ref have 0 counts for %s and %s", gene_id, var_id)
                continue
            total_coverage = hap_count_subset['refCount'].sum() + hap_count_subset['altCount'].sum()
            real_counts = [[hap_count_subset['refCount']], [hap_count_subset['altCount']]]
            # chi2, p_value_real, dof, expected = chi2_contingency(real_counts)
            chi2, p_value_real, _, _ = chi2_contingency(real_counts)
            cohen = round(numpy.sqrt(chi2 / number_of_isoforms), 8)
            data.append([gene_id, var_id, total_coverage, number_of_isoforms, cohen, p_value_real])
            output_power = pandas.DataFrame({
                'gene_id': [item[0] for item in data], 'variant_id': [item[1] for item in data],
                'total_coverage': [item[2] for item in data], 'number_of_isoforms': [item[3] for item in data],
                'cohen': [item[4] for item in data], 'p_value': [item[5] for item in data]
            })
            # select the SNP with the highest sum count in case of multiple SNP for a gene
            output_power_filtered = (output_power.groupby("gene_id").apply(_f))
            output_power_filtered_indexed = output_power_filtered.set_index('gene_id')
            # output_power_filtered_indexed['qvalue'] = fdr_adjust(output_power_filtered_indexed['p_value'].tolist())
            output_power_filtered_indexed['qvalue'] = pvalue_adjust(pvalues=output_power_filtered_indexed['p_value'].tolist())
            return output_power_filtered_indexed
            # output_power_filtered_indexed.to_csv(output_dir + "processed_asts_" + file, sep='\t')

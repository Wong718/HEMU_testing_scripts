#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project ：HEMU_Database_Main
@File    ：tasks.py
@Author  ：Edward Zhu, zhuyzh37@mail2.sysu.edu.cn
@Date    ：2022/10/2 18:00
@IDE     ：PyCharm
-----------------------------
Description: Celery-dependent async tasks
Input:
Output:

Environment requirement (if any): HEMUdb_new (Remote Server, Ubuntu 20.04)
- Dependencies: Django(v3.1), Project files, etc.
-----------------------------
Notes:

-----------------------------
File Revisions:
    2022/10/2: Version 1 - Creation
    2023/2/18: Version 2 - Finished async conversion for all transcriptome-derived analysis
"""

# ===================================================================================
# Global imports
# Async-dependent
from __future__ import absolute_import, unicode_literals
from celery import shared_task
# Function-dependent
import time
import pandas as pd
import random
import os
# ===================================================================================

# ===================================================================================
# Module: Transcriptome-derived analysis
# Gene Expression Query
from Mainapp.Main_scripts.transcriptome import GeneExpressionDataCurator
from Mainapp.R_visualization.transcriptome import gene_expression_query_plt_generator
# Differential Gene Expression Analysis
from Mainapp.Main_scripts.transcriptome import GeneDEDataCurator
from Mainapp.R_visualization.transcriptome import gene_differential_analysis_plt_generator
# GO / KEGG Enrichment
from Mainapp.Main_scripts.transcriptome import GoKeggDataCurator
from Mainapp.R_visualization.transcriptome import gene_gokegg_enrichment_plt_generator
# WGCNA
from Mainapp.Main_scripts.transcriptome import WGCNASampleGeneDataCurator
from Mainapp.R_visualization.transcriptome import gene_wgcna_section_deployer
# Sequence Obtain
from Mainapp.Main_scripts.transcriptome import GeneSequenceObtainer
# ===================================================================================

# ===================================================================================
# Module: Transcriptome-derived analysis
# Phylogenetic Analysis
from Mainapp.Main_scripts.genefamily import PhyloAnalysisHandler
from Mainapp.R_visualization.genefamily import genefam_phylogenetic_plt_generator
# Family Identification - hmmsearch & blastp
from Mainapp.Main_scripts.genefamily import GeneFamSearchHandler

# ===================================================================================
# Module: Multi-omics anaysis
# miRNA-seq
from Mainapp.Main_scripts.miRNA import miRNAExpressionDataCurator
from Mainapp.R_visualization.miRNA import miRNA_expression_query_plt_generator

# ChIP-seq
from Mainapp.R_visualization.ChIPseq import ChIP_analysis_plt_generator

# Sample task to validate if celery backend is working
@shared_task
def testing_task1():
    time.sleep(10)
    pass


# ===================================================================================
# Module: Transcriptome-derived analysis
@shared_task
def gene_expression_query_handler(query_species_name, exp_sheet_name, gene_list, expression_format):
    """

    :param query_species_name:
    :param exp_sheet_name:
    :param gene_list:
    :param expression_format:
    :param var_transfer_key:
    :return:
    """
    query_list_full = []  # [[gene_ID, expressed_samples, total_samples, max, min, median], [..]]
    for indv_gene in gene_list:
        try:
            # Curate expression dataframe
            exp_df = GeneExpressionDataCurator.gene_exp_df_builder(indv_gene, exp_sheet_name)
            # Add basic expression statistics (expressed samples, avg, median, etc.)
            query_list_full.append(
                GeneExpressionDataCurator.fund_info_obtainer(indv_gene, exp_df, expression_format)
            )
            # Generate plots
            gene_expression_query_plt_generator.overview_barplot(
                exp_df, indv_gene, expression_format, query_species_name
            )
            gene_expression_query_plt_generator.tissue_specific_barplot(
                exp_df, indv_gene, expression_format, query_species_name
            )
            # Generate gene-specific expression table (CSV format)
            exp_df.to_csv("Mainapp/static/Temp_R_html/" + indv_gene + "_values.csv")

        except:
            return RuntimeError
    # Returning full gene expression information list to views function
    # [query_list, query_list_full, query_format, query_species]
    gene_exp_query_return_list = [gene_list, query_list_full, expression_format, query_species_name]
    return gene_exp_query_return_list


@shared_task
def gene_differential_analysis_handler(query_species_name, exp_sheet_name,
                                       group1_samples_list, group1_name, group2_samples_list, group2_name,
                                       logfc_threshold, pvalue_threshold, heatmap_gene_count):
    """

    :param query_species_name:
    :param exp_sheet_name:
    :param group1_samples_list:
    :param group1_name: [list]
    :param group2_samples_list:
    :param group2_name: [list]
    :param logfc_threshold:
    :param pvalue_threshold:
    :param heatmap_gene_count:
    :return:
    """
    try:
        DE_data_raw, DE_group_list, DE_group_color_list = GeneDEDataCurator.gene_de_df_builder(
            exp_sheet_name, group1_samples_list, group1_name, group2_samples_list, group2_name
        )
        task_destination_folder = gene_differential_analysis_plt_generator.GeneDifferentialAnalysis(
            DE_data_raw, DE_group_list, DE_group_color_list,
            logfc_threshold, pvalue_threshold, heatmap_gene_count,
            group1_name, group2_name
        )
    except:
        return RuntimeError
    # Returning full differential gene expression information list to views function
    gene_differential_analysis_return_list = [task_destination_folder, query_species_name,
                                              group1_name, group1_samples_list, group2_name, group2_samples_list,
                                              logfc_threshold, pvalue_threshold, heatmap_gene_count]
    return gene_differential_analysis_return_list


@shared_task
def gene_gokegg_enrichment_handler(gene_list, gokegg_sheet_name, enrich_format, species_filename):
    try:
        gokegg_df = GoKeggDataCurator.gene_gokegg_df_builder(gene_list, gokegg_sheet_name)

        plt_filename_id = gene_gokegg_enrichment_plt_generator.gokegg_enrich_plt(
            gokegg_df, enrich_format, species_filename
        )
    except:
        return RuntimeError
    gene_gokegg_enrichment_return_list = [plt_filename_id]
    return gene_gokegg_enrichment_return_list


@shared_task
def gene_wgcna_handler(args_dict, section_number):
    """

    :param args_dict:
    :param section_number:
    :return:
    """
    #try:
    if section_number == "1":

        # Generate unique folder name and create the folder
        dirname = "wgcna" + str(random.randint(int(1e8), int(1e9) - 1))
        os.makedirs("Mainapp/static/Temp_R_wgcna/" + dirname, exist_ok=True)

        # exp_sheet_name, gene_list, sample_list, expression_format
        sample_exp_df = WGCNASampleGeneDataCurator.wgcna_init_df_builder(
            args_dict["query_table"],
            args_dict["gene_list_final"],
            args_dict["accession_list_final"],
            args_dict["query_format"],
        )
        sampleGeneData = "Mainapp/static/Temp_R_wgcna/" + dirname + "/sample_expression_data.csv"
        sample_exp_df.to_csv(sampleGeneData)

        # sampleGeneData, RcCutoff, samplePerc, GeneNum, cutmethod, rscut, datatype, anamethod, dirname
        dirname = gene_wgcna_section_deployer.wgcna_section1_deployer(
            sampleGeneData,
            args_dict["exp_cut_threshold"],
            args_dict["sample_cut_threshold"],
            args_dict["genes_retain"],
            args_dict["cut_method"],
            args_dict["soft_power_cutoff"],
            args_dict["query_format"],
            args_dict["transformation_method"],
            dirname,
        )

        return [dirname]
    #except:
    #    return RuntimeError


@shared_task
def gene_sequence_obtain_handler(query_species_name_filename, gene_list, query_species_name):
    try:
        # return a list containing all the fasta entries, with each element representing a single fasta line.
        sequence_fasta_filename = GeneSequenceObtainer.gene_sequence_query(gene_list, query_species_name_filename,
                                                                           query_species_name, True)
    except:
        return RuntimeError
    gene_sequence_obtain_return_list = [sequence_fasta_filename]
    return gene_sequence_obtain_return_list


# ===================================================================================

# ===================================================================================
# Gene family analysis
@shared_task
def genefam_phylogenetic_analysis_handler(query_species_name, query_sequences,
                                          msa_method, pairwise_dist_method, tree_layout_method, bootstrap_rep_num):
    try:
        sequence_fasta_string = PhyloAnalysisHandler.genefam_data_validation(query_sequences, query_species_name)
        result_folder_name = genefam_phylogenetic_plt_generator.genefam_phylogenetic_analysis(
            sequence_fasta_string,
            msa_method,
            pairwise_dist_method,
            tree_layout_method,
            bootstrap_rep_num
        )
    except:
        return RuntimeError
    genefam_phylo_analysis_return_list = [result_folder_name]
    return genefam_phylo_analysis_return_list


@shared_task
def genefam_identification_hmmsearch_handler(species_filename, hmm_query_filename,
                                             seq_evalue_threshold, dom_evalue_threshold):
    try:
        result_folder_name = GeneFamSearchHandler.hmmsearch_deployer(species_filename, hmm_query_filename,
                                                                     seq_evalue_threshold, dom_evalue_threshold)
        return [result_folder_name]
    except:
        return RuntimeError


# ===================================================================================
# ===================================================================================
# Module: Multiomics data analysis
# miRNA-seq
@shared_task
def miRNA_expression_query_handler(query_species_name, exp_sheet_name, miRNA_list, expression_format):
    """

    :param query_species_name:
    :param exp_sheet_name:
    :param miRNA_list:
    :param expression_format:
    :param var_transfer_key:
    :return:
    """
    query_list_full = []  # [[miRNA_ID, expressed_samples, total_samples, max, min, median], [..]]
    for indv_gene in miRNA_list:
        try:
            # Curate expression dataframe
            exp_df = miRNAExpressionDataCurator.gene_exp_df_builder(indv_gene, exp_sheet_name)
            # Add basic expression statistics (expressed samples, avg, median, etc.)
            query_list_full.append(
                miRNAExpressionDataCurator.fund_info_obtainer(indv_gene, exp_df, expression_format)
            )
            # Generate plots
            miRNA_expression_query_plt_generator.overview_barplot(
                exp_df, indv_gene, expression_format, query_species_name
            )
            miRNA_expression_query_plt_generator.tissue_specific_barplot(
                exp_df, indv_gene, expression_format, query_species_name
            )
            # Generate gene-specific expression table (CSV format)
            exp_df.to_csv("Mainapp/static/Temp_R_miRNA/" + indv_gene + "_values.csv")

        except:
            return RuntimeError
    # Returning full gene expression information list to views function
    # [query_list, query_list_full, query_format, query_species]
    miRNA_exp_query_return_list = [miRNA_list, query_list_full, expression_format, query_species_name]
    return miRNA_exp_query_return_list

# ChIP-seq
@shared_task
def ChIP_query_handler(sample_id, tssRegion, species,ignore_1st_exon, ignore_1st_intron, ignore_downstream,ignore_promoter_subcategory):
    """

    :param sample_id:
    :param tssRegion:
    :param species:
    :param ignore_1st_exon:
    :param ignore_1st_intron:
    :param ignore_downstream:
    :param ignore_promoter_subcategory:
    :return:
    """
    try:
        # Generate plots and files
        output_folder_name = ChIP_analysis_plt_generator.SingleChIPAnalysis(
            sample_id, tssRegion, species,ignore_1st_exon = 'F', ignore_1st_intron = 'F', ignore_downstream = 'F',ignore_promoter_subcategory='F'
        )

    except:
        return RuntimeError

    ChIP_analysis_return_list = [sample_id, tssRegion, species, ignore_1st_exon, ignore_1st_intron, ignore_downstream, ignore_promoter_subcategory,output_folder_name]
    return ChIP_analysis_return_list
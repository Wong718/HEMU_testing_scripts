import pandas as pd
import numpy as np
import pymysql
from Mainapp.Main_scripts.MainConfiguration import query_sql

# Extract gene expression data from database and append tissue information regarding each sample
# Database configurations

# Used for pymysql services
dbhost = str(query_sql("host"))
dbuser = str(query_sql("user"))
dbpassword = str(query_sql("pwd"))
dbdatabase = str(query_sql("dbname"))


def gene_exp_df_builder(miRNA_id, exp_sheet_name):
    """

    :param miRNA_id:
    :param exp_sheet_name:
    :param sampleinfo_sheet_name:
    :return: pd.DataFrame object, with column index: [sample_id, Counts, Normalized, sample_tissue]
    """
    global colnames, results, init_df_list, result_tissue

    init_df_list = []  # List for generating initial dataframe (sample_id , Counts, Normalized)
    db = pymysql.connect(host=dbhost, user=dbuser, password=dbpassword, database=dbdatabase)
    cursor = db.cursor()

    sqlcmd_select_gene = "SELECT * FROM %s WHERE gene='%s';" % (exp_sheet_name, miRNA_id)

    try:  # Execute SQL command
        cursor.execute(sqlcmd_select_gene)
        results = cursor.fetchall()
    except:
        print("Exception occurred while querying database.")

    # Build dataframe
    init_df_list = [[indv[1], indv[2], indv[3]] for indv in results]
    init_df = pd.DataFrame(init_df_list, columns=['sample_id', 'Counts', 'Normalized', 'tissue_type'])
    # SRR-ID COUNTS NORMALIZED TISSUE INFO1 INFO2 INFO3

    return init_df


def fund_info_obtainer(miRNA_id, gene_exp_df, exp_format):
    """
    Obtain critical information regarding gene expression, Used in the right-side panel of miRNA-seq search result.
    :param miRNA_id: miRNA ID
    :param gene_exp_df: sample_id Counts Normalized tissue_type
    :param exp_format: Counts/Normalized, used for extracting expression level from dataframe
    :return: gene_exp_detail = [miRNA_id, expressed_samples, total_samples, max, min, median]
    """
    gene_exp_detail = []
    expressed_threshold = 1  # Modifiable, default with both Counts and Normalized.

    if exp_format == "Counts":
        exp_level_comp = gene_exp_df['Counts'].tolist()
        gene_exp_detail.append([
            miRNA_id,
            len([indv for indv in exp_level_comp if indv >= expressed_threshold]),
            len(exp_level_comp),
            max(exp_level_comp),
            min(exp_level_comp),
            np.median(exp_level_comp),
        ])
    elif exp_format == "Normalized":
        exp_level_comp = gene_exp_df['Normalized'].tolist()
        gene_exp_detail.append([
            miRNA_id,
            len([indv for indv in exp_level_comp if indv >= expressed_threshold]),
            len(exp_level_comp),
            max(exp_level_comp),
            min(exp_level_comp),
            np.median(exp_level_comp),
        ])

    return gene_exp_detail

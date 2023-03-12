import pandas as pd
import numpy as np
import pymysql
from Mainapp.Main_scripts.MainConfiguration import query_sql, query_tables
from Mainapp.Main_scripts.TFHeatmapSampleConfiguration import return_sample

# Used for pymysql services
dbhost = str(query_sql("host"))
dbuser = str(query_sql("user"))
dbpassword = str(query_sql("pwd"))
dbdatabase = str(query_sql("dbname"))


def TF_fam_to_geneid_query(species, TF_fam_name_list, indv_geneid_list):
    """
    :param TF_fam_name_list: single gene id, not accompanied by a family (list)
    :param TF_fam_name_list: [YABBY, MYB]
    :return: _TF_fam_list: [[TF_geneid, TF_family], [TF_geneid, TF_family], ...]
    """
    global results
    _TF_fam_list = []
    db = pymysql.connect(host=dbhost, user=dbuser, password=dbpassword, database=dbdatabase)
    cursor = db.cursor()

    for indv_TF_fam in TF_fam_name_list:
        sqlcmd_select_tf= "SELECT * FROM %s WHERE family='%s';" % \
                                      (query_tables(species + "_tf"), indv_TF_fam)
        try:  # Execute SQL command
            cursor.execute(sqlcmd_select_tf)
            results = cursor.fetchall()
        except:
            print("Exception occurred.")
        for indv_index in results:
            _TF_fam_list.append([indv_index[1], indv_index[2]])  # [TF_geneid, TF_family]

    for indv_TF_gene in indv_geneid_list:
        _TF_fam_list.append([indv_TF_gene, 'custom_TF'])

    return _TF_fam_list  # Embedded list, not a dataframe


def TF_heatmap_df_builder(species, TF_fam_list, tissue_list):
    """

    :param species:
    :param TF_fam_list: [[TF_geneid, TF_family], [TF_geneid, TF_family], ...]
    :param tissue_list: [leaf, root], to be matched with query against samples
    :return: dataframe, [sample_id, sample_tissue, gene_id, tf_fam, tpm]
    """

    main_TF_heatmap_res = [] # [sample_id, sample_tissue, gene_id, tf_fam, tpm]
    db = pymysql.connect(host=dbhost, user=dbuser, password=dbpassword, database=dbdatabase)
    cursor = db.cursor()

    for indv_tissue in tissue_list:
        # print("cataloging %s" % indv_tissue, return_sample(species)[indv_tissue])
        for indv_sample in return_sample(species)[indv_tissue]:  # list of SRR samples
            for indv_tf in TF_fam_list:
                # print(indv_tf)
                sqlcmd_select_tf_expression = "SELECT * FROM %s WHERE gene='%s' and sample_id='%s';" % \
                                              (query_tables(species + "_exp"), indv_tf[0], indv_sample)
                # print(sqlcmd_select_tf_expression)
                try:  # Execute SQL command
                    cursor.execute(sqlcmd_select_tf_expression)
                    sql_result_tpm = str(cursor.fetchall()[0][3])
                    # print(sql_result_tpm)
                    main_TF_heatmap_res.append([indv_sample,
                                                indv_tissue,
                                                indv_tf[0],
                                                indv_tf[1],
                                                str(sql_result_tpm)
                                                ])
                except:
                    print("Exception occurred while querying database.")

    # Build dataframe
    init_df = pd.DataFrame(main_TF_heatmap_res,
                           columns=["sample_id", "sample_tissue", "gene_id", "tf_fam", "tpm"])
    return init_df

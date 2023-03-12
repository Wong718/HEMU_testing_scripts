"""
MainConfiguration.py
Stores critical information and table names used in the project
"""
tables_dict = {  # Table names
    'coix_exp': 'coix_exp',
    'coix_samp': 'coix_samp',
    'coix_tf': 'coix_tf',
    'coix_te': 'coix_te',
    'coix_gene_raw': 'coix_gene_raw',
    'coix_gokegg': 'coix_gokegg',
    'coix_filename': 'Clacr',

    'zea_exp': 'zea_exp',
    'zea_samp': 'zea_samp',
    'zea_tf': 'zea_tf',
    'zea_te': 'zea_te',
    'zea_gene_raw': 'zea_gene_raw',
    'zea_gokegg': 'zea_gokegg',
    'zea_filename': 'Zmaysv4',

    'sorghum_exp': 'sorghum_exp',
    'sorghum_samp': 'sorghum_samp',
    'sorghum_tf': 'sorghum_tf',
    'sorghum_gene_raw': 'sorghum_gene_raw',
    'sorghum_gokegg': 'sorghum_gokegg',
    'sorghum_filename': 'Sbico',

    'saccharum_exp': 'saccharum_exp',
    'saccharum_samp': 'saccharum_samp',
    'saccharum_gene_raw': 'saccharum_gene_raw',
    'saccharum_filename': 'Sspon',

    'miscanthus_filename': 'Mluta',

    'chrysopogon_filename': 'Cseru',

    'hyparrhenia_filename': 'Hdipla',

    'themeda_filename': 'Ttria',
}

sql_user = ""
sql_pwd = ""
sql_host = "localhost"
sql_port = ""
sql_dbname = ""

admin_pwd = ""


def query_tables(tbl2query):
    try:
        return tables_dict[tbl2query]
    except NameError:
        return None


def query_admin_pwd():
    return admin_pwd


def query_sql(identity):
    if identity == "host":
        return sql_host
    elif identity == "user":
        return sql_user
    elif identity == "pwd":
        return sql_pwd
    elif identity == "port":
        return sql_port
    elif identity == "dbname":
        return sql_dbname
    else:
        return None

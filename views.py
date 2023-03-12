#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project ：HEMU_Database_Main
@File    ：views.py
@Author  ：Edward Zhu, zhuyzh37@mail2.sysu.edu.cn
@Date    ：2022/8/27 18:00
@IDE     ：PyCharm
-----------------------------
Description: The main views function for HEMU database.
Input:
Output:

Environment requirement (if any): HEMUdb_new (Remote Server, Ubuntu 20.04)
- Dependencies: Django(v3.1), Project files, etc.
-----------------------------
Notes:

-----------------------------
File Revisions:
    2022/8/27: Version 1 - Creation
"""
import json

from django.shortcuts import render, HttpResponse, redirect

# Server Configurations
from Mainapp.Main_scripts import MainConfiguration
from Mainapp.Main_scripts.site_manage import TmpFileCleaner

# TF Module
from Mainapp.Main_scripts import TFHeatmapDataCurator
from Mainapp.R_visualization import TF_heatmap_generator
from Mainapp.Main_scripts import TFHeatmapSampleConfiguration

# TE Module
from Mainapp.Main_scripts import TEBySampleDataCurator
from Mainapp.R_visualization import TE_expression_plt_generator
from Mainapp.Main_scripts import TEByFamilyDataCurator

# Jbrowse2 Module
from Mainapp.Main_scripts.jbrowse import jbrowse_species_catalog

# Celery applications
import Mainapp.tasks as celery_tasks  # Async task list
from celery.result import AsyncResult  # Interface for task status querying


def get_async_task_progress(request):
    """
    Global view function for querying task progeess.
    :param request: AJAX request, GET, ?taskid=[celery_task_id]
    :return: json, {'state': 'PENDING'} / {'state': 'SUCCESS'} / {'state': 'FAILURE'}
    """
    if request.is_ajax():
        taskid = request.GET.get('taskid')
        if taskid:
            result = AsyncResult(taskid)
            # Another validation method, print(result.successful())
            response_data = {'state': result.state, }
            return HttpResponse(json.dumps(response_data), content_type='application/json')
    return render(request, 'error.html', {'error_message': 'INVALID_REQUEST'})


def render_asyc_onhold(request):
    task_id = request.GET.get('taskid')
    redirect_url = request.GET.get('dir')
    if task_id and redirect_url:
        return render(request, 'async_task_onhold.html',
                      {'task_id': task_id, 'redirect_url': redirect_url})
    return render(request, 'error.html', {'error_message': 'INVALID_URL_QUERY'})


def async_test(request):
    """
    Testing task deployer, submit a 10-second time consuming task to celery queue.
    :param request:
    :return:
    """
    if request.method == "GET":
        if request.GET.get('task') == "task1":
            async_result = celery_tasks.testing_task1.delay()
            return render(request, 'async_test/test_progressbar.html', {'task_id': async_result.task_id})
        elif request.GET.get('success') == "1":
            return render(request, 'error.html', {'error_message': 'TASK_IS_COMPLETED'})

        return render(request, 'async_test/test_progressbar.html')


# ==============MODULE: HOME PAGE==============

# 1.2 HOME PAGE RENDERER
def init_scr(request):
    return render(request, 'home.html')


# 1.3 SITE MANAGER VALIDATOR & RENDERER
def site_manage(request):
    if request.method == "GET":
        # Common GET, request main page
        return render(request, 'site_manage/site_manage_auth.html')
    elif request.method == "POST":

        if request.POST.get('clear_tmp_files') == "1":
            identity_number_deleted = int(TmpFileCleaner.clean_tmp_files())
            tmp_message = "Cleaning completed, removed %d identities." % identity_number_deleted
            return render(request, 'site_manage/site_manage_dashboard.html',
                          {'clear_tmp_files_message': tmp_message})

        # Authentication
        if str(request.POST.get('password')) == MainConfiguration.query_admin_pwd():
            return render(request, 'site_manage/site_manage_dashboard.html')
        else:
            return render(request, 'site_manage/site_manage_auth.html',
                          {'error_message': 'Password Incorrect!'})


# ==============MODULE: TRANSCRIPTOME-DERIVED ANALYSIS==============
# 2.1 GENE EXPRESSION PROFILE QUERY HANDLER
def gene_expression_query_async(request):
    if request.method == 'GET':
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Task completed, return [query_list, query_list_full, query_format, query_species]
                    context = {
                        'query_list': task_return[0],
                        'query_list_full': task_return[1],
                        'query_format': task_return[2],
                        'query_species': task_return[3],
                        'last_query': ";".join(task_return[0]),
                    }
                    return render(request, 'gene_expression/gene_search_main.html', context)

                # task_return == None means no database queue hit, probably due to faulty input.
                return render(request, 'gene_expression/gene_search_main.html',
                              {'error_message': 'Invalid entry. Please check gene nomenclature and submit the query '
                                                'again.'})

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE'})
        # Render default view
        return render(request, 'gene_expression/gene_search_main.html')

    elif request.method == "POST":
        # Task submission
        # Obtain query from frontend form submission
        query_list = request.POST.get('main_query').split(';')  # Original query data
        query_species = request.POST.get('query_species')  # Species to query
        query_format = request.POST.get('optionsRadios')  # FPKM / TPM

        # Obtain sql tables
        query_exp_table = MainConfiguration.query_tables(str(query_species) + "_exp")

        # Deploy task
        async_result = celery_tasks.gene_expression_query_handler.delay(
            query_species,
            query_exp_table,
            query_list,
            query_format,
        )
        # Return page containing task id, awaiting frontend recognition
        return render(request, 'gene_expression/gene_search_main.html', {'task_id': async_result.task_id})


# 2.2 SEQUENCE ACQUISITION: Gene, Transcript, Protein
def gene_sequence_obtain_async(request):
    if request.method == "GET":

        # Submit GET-format task, used for small amount of querying.
        query_species = request.GET.get("sp")
        query_raw = request.GET.get("gene")
        query_format = request.GET.get("format")
        if query_species and query_raw and query_format:
            if len(query_raw.split(";")) > 10:
                return render(request, 'error.html', {'error_message': 'QUERY_TOO_LONG'})

            # Obtain filenames for fasta file containing gene sequences
            species_filename = MainConfiguration.query_tables(query_species + "_filename")
            # Deploy task
            async_result = celery_tasks.gene_sequence_obtain_handler.delay(
                species_filename, query_raw.split(";"), query_format
            )
            # Return async task manager, awaiting frontend recognition
            return render(request, 'gene_expression/gene_sequence_acquire.html', {'task_id': async_result.task_id})

        # Process task return status
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Task completed, return list
                    # [sequence_fasta_filename]
                    if not task_return[0].startswith("sequence"):
                        # Empty result list
                        # Nomenclature: "sequence" + str(random.randint(int(1e8), int(1e9) - 1)) + ".fasta"
                        return render(request, 'gene_expression/gene_sequence_acquire.html',
                                      {'error_message': 'Invalid entry. Please check gene nomenclature and submit the '
                                                        'query again.'})

                    # Read file content and render to frontend
                    try:
                        with open("Mainapp/static/Temp_R_genefam/" + str(task_return[0]), mode='r') as fasta_fh:
                            fasta_content = fasta_fh.read()

                    except FileNotFoundError:
                        return render(request, 'gene_expression/gene_sequence_acquire.html',
                                      {'error_message': 'Invalid entry. Please check gene nomenclature and submit the '
                                                        'query again.'})
                    context = {
                        'result_list_full': [fasta_content],
                        'fasta_filename': task_return[0],
                    }
                    return render(request, 'gene_expression/gene_sequence_acquire.html', context)
                else:
                    return render(request, 'error.html', {'error_message': 'EMPTY_TASK_RESPONSE'})

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE'})
        # Render default view
        return render(request, 'gene_expression/gene_sequence_acquire.html')

    elif request.method == "POST":

        # Redirect for BLAST instance
        seq_transfer = request.POST.get("sequence_raw")
        if seq_transfer:
            return render(request, 'BLAST/sequenceserver_display.html', {'seq_transfer': seq_transfer, })

        # Task submission
        # Obtain query from frontend form submission
        query_list = request.POST.get('main_query').split(';')  # Original query data
        query_species = request.POST.get('query_species')  # Species to query
        query_format = request.POST.get('optionsRadios')  # gene / transcript / protein

        # Obtain filenames for fasta file containing gene sequences
        species_filename = MainConfiguration.query_tables(query_species + "_filename")
        # Deploy task (write .fasta file to disk in order to prevent celery overflow)
        async_result = celery_tasks.gene_sequence_obtain_handler.delay(
            species_filename, query_list, query_format
        )
        # Return async task manager, awaiting frontend recognition
        return render(request, 'gene_expression/gene_sequence_acquire.html', {'task_id': async_result.task_id})


# 2.3 DGE QUERY HANDLER
def gene_differential_analysis_async(request):
    if request.method == "GET":
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Task completed, return list
                    # [task_destination_folder, query_species_name, group1_name, group1_samples_list,
                    # group2_name, group2_samples_list, logfc_threshold, pvalue_threshold, heatmap_gene_count]
                    context = {
                        'task_destination_folder': task_return[0],
                        'species_query': task_return[1],
                        'group1_name': task_return[2],
                        'group1_samples_list': task_return[3],
                        'group2_name': task_return[4],
                        'group2_samples_list': task_return[5],
                        'logfc_threshold': task_return[6],
                        'pvalue_threshold': task_return[7],
                        'heatmap_gene_count': task_return[8],
                    }
                    return render(request, 'gene_expression/gene_DE_report_display.html', context)

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE_EXCEPTION'})
        # Render default view
        return render(request, 'gene_expression/gene_DE_main.html')

    elif request.method == "POST":
        # Task submission
        # Obtain query from frontend form submission
        query_species = request.POST.get('species_query')
        logfc_threshold = request.POST.get('logfc_threshold')
        pvalue_threshold = request.POST.get('pvalue_threshold')
        heatmap_gene_count = request.POST.get('heatmap_gene_count')

        group1_samples_list = request.POST.get('group1_samples').split(";")
        group1_name = request.POST.get('group1_name')
        group2_samples_list = request.POST.get('group2_samples').split(";")
        group2_name = request.POST.get('group2_name')

        # Obtain sql tables
        query_exp_table = MainConfiguration.query_tables(str(query_species) + "_exp")

        # Deploy task
        async_result = celery_tasks.gene_differential_analysis_handler.delay(
            query_species, query_exp_table,
            group1_samples_list, group1_name, group2_samples_list, group2_name,
            logfc_threshold, pvalue_threshold, heatmap_gene_count
        )
        # Return async task manager, awaiting frontend recognition
        return redirect("/HEMUdb/asynctasks/onhold?taskid=%s&dir=%s" %
                        (async_result.task_id, "/HEMUdb/gene/dge"))


# 2.3 DGE STATIC FILE LOAD HANDLER
def load_DE_staticfile(request, identifier_name, file_name, frame_height, frame_width):
    return render(request, 'static_html_display.html',
                  {
                      'identifier_name': identifier_name,
                      'file_name': file_name,
                      'frame_height': frame_height,
                      'frame_width': frame_width,
                  })


# 2.4 GO/KEGG ENRICHMENT
def gene_gokegg_enrichment_async(request):
    if request.method == "GET":
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Task completed, return list
                    # [plt_filename_id]
                    context = {
                        'filename_id': task_return[0],
                    }
                    return render(request, 'gene_expression/gokegg_enrich_result.html', context)

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE_EXCEPTION'})
        # Render default view
        return render(request, 'gene_expression/gokegg_enrich.html')

    elif request.method == "POST":
        # Task submission
        # Obtain query from frontend form submission
        species_query = request.POST.get("query_species")
        enrich_format = request.POST.get('optionsRadios')  # GO / KEGG
        # Separate query and remove blank entries
        gene_list_final = [indv_gene.rstrip("\r") for indv_gene in request.POST.get("query_gene_list").split("\n")]
        gene_list_final = [indv_entry for indv_entry in gene_list_final if indv_entry != '']

        # Obtain sql tables
        query_table = MainConfiguration.query_tables(str(species_query) + "_gokegg")
        # Obtain filename for accessing background annotation
        species_filename = MainConfiguration.query_tables(str(species_query) + "_filename")

        # Deploy task
        # params: gene_list, gokegg_sheet_name, enrich_format
        async_result = celery_tasks.gene_gokegg_enrichment_handler.delay(
            gene_list_final, query_table, enrich_format, species_filename
        )
        # Return async task manager, awaiting frontend recognition
        return redirect("/HEMUdb/asynctasks/onhold?taskid=%s&dir=%s" %
                        (async_result.task_id, "/HEMUdb/gene/gokegg"))


# 2.5 WGCNA ANALYSIS
def gene_wgcna_shinyappgen_async(request):
    if request.method == "GET":
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Task completed, return list
                    # [output_app_dir, output_app_dir_name]
                    # context = {
                    #    'shiny_app_name': task_return[1],
                    # }
                    # return render(request, 'gene_expression/wgcna_app.html', context)
                    return HttpResponse("Task completed, %s" % str(task_return))

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE_EXCEPTION'})
        # Render default view (data submission)
        return render(request, 'gene_expression/wgcna_ana_dashboard.html')

    elif request.method == "POST":
        # Task submission
        # Create a new instance, obtain query from frontend form submission
        form_section_id = request.POST.get("form_section_id")

        if form_section_id == "1":
            query_species = request.POST.get("query_species")
            query_format = request.POST.get("query_exp_met")  # FPKM / TPM
            transformation_method = request.POST.get("transformation_method")  # rawFPKM / logFPKM
            exp_cut_threshold = request.POST.get("exp_cut_threshold")
            sample_cut_threshold = request.POST.get("sample_cut_threshold")
            cut_method = request.POST.get("cut_method")
            genes_retain = request.POST.get("genes_retain")
            soft_power_cutoff = request.POST.get("soft_power_cutoff")

            # Obtain gene list
            gene_list_final = [indv_gene.rstrip("\r") for indv_gene in request.POST.get("query_gene_list").split("\n")]
            gene_list_final = [indv_entry for indv_entry in gene_list_final if indv_entry != '']  # Remove blank entries
            # Obtain accession list
            accession_list_final = [indv_accession.rstrip("\r") for indv_accession in
                                    request.POST.get("query_accession_list").split("\n")]
            accession_list_final = [indv_entry for indv_entry in
                                    accession_list_final if indv_entry != '']  # Remove blank entries
            # Check gene list & accession list length
            if not 1 < len(gene_list_final) < 4001:
                return render(request, 'gene_expression/wgcna_ana_dashboard.html', {
                    'error_message': 'Invalid gene list length, please re-submit the query.'
                })
            if not 1 < len(accession_list_final) < 61:
                return render(request, 'gene_expression/wgcna_ana_dashboard.html', {
                    'error_message': 'Invalid sample accession list length, please re-submit the query.'
                })

            # Obtain sql tables for expression data
            query_table = MainConfiguration.query_tables(str(query_species) + "_exp")
            # Pack function context
            task_context = {
                'query_table': query_table,
                'query_species': query_species,
                'query_format': query_format,
                'transformation_method': transformation_method,
                'exp_cut_threshold': exp_cut_threshold,
                'sample_cut_threshold': sample_cut_threshold,
                'cut_method': cut_method,
                'genes_retain': genes_retain,
                'soft_power_cutoff': soft_power_cutoff,
                'gene_list_final': gene_list_final,
                'accession_list_final': accession_list_final,
            }
            async_result = celery_tasks.gene_wgcna_handler.delay(task_context, "1")
            return redirect("/HEMUdb/asynctasks/onhold?taskid=%s&dir=%s" %
                            (async_result.task_id, "/HEMUdb/gene/wgcna?section=1\&"))

        elif form_section_id == "2":
            query_project_id = request.POST.get("query_projectid")
            query_sftpower = request.POST.get("query_sftpower")

            return HttpResponse("Sub section 2, %s %s" % (query_project_id, query_sftpower))

        elif form_section_id == "3":
            query_project_id = request.POST.get("query_projectid")
            sampletrait_tbl_file = request.FILES.get("query_sampletrait_tbl_file", None)
            return HttpResponse("Sub section 3, %s" % query_project_id)

        # async_result = celery_tasks.gene_WGCNA_handler.delay(
        #    query_table,
        #    gene_list_final,
        #    accession_list_final,
        #    query_format,
        # )
        # Return async task manager, awaiting frontend recognition
        # return redirect("/HEMUdb/asynctasks/onhold?taskid=%s&dir=%s" %
        #                (async_result.task_id, "/HEMUdb/gene/wgcna"))


# ============MODULE: GENE FAMILY ANALYSIS=============

def genefam_phylogenetic_analysis_async(request):
    if request.method == "GET":
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Task completed, return list
                    # [plt_filename_id]
                    context = {
                        'result_folder_name': task_return[0],
                    }
                    return render(request, 'gene_family/gene_phylogenic_analysis_result.html', context)

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE_EXCEPTION'})
        # Render default view
        return render(request, "gene_family/gene_phylogenetic_analysis.html")

    elif request.method == "POST":
        # Task submission
        # Obtain query from frontend form submission
        species_query = request.POST.get("query_species")  # by HTML tag name
        sequence_query = request.POST.get("query_sequence")
        msa_method = request.POST.get("query_msa_method")
        pairwise_dist_method = request.POST.get("query_pairwise_dist_method")
        tree_layout_method = request.POST.get("query_tree_layout_method")
        bootstrap_rep_num = request.POST.get("query_bootstrap_rep_num")
        try:
            # Validate if variable bootstrap_rep_num is a number
            int(bootstrap_rep_num) / 2
        except ValueError:
            return render(request, 'error.html', {'error_message': 'INVALID_QUERY'})

        # Deploy task
        # params: gene_list, gokegg_sheet_name, enrich_format
        async_result = celery_tasks.genefam_phylogenetic_analysis_handler.delay(
            species_query, sequence_query,
            msa_method, pairwise_dist_method, tree_layout_method, bootstrap_rep_num
        )
        # Return async task manager, awaiting frontend recognition
        return redirect("/HEMUdb/asynctasks/onhold?taskid=%s&dir=%s" %
                        (async_result.task_id, "/HEMUdb/genefam/phylo"))


def genefam_identification_hmm_async(request):
    if request.method == "GET":
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Task completed, return list
                    # [plt_filename_id]
                    if task_return[0] == 102:
                        return render(request, 'gene_family/gene_hmmsearch_analysis.html',
                                      {'error_message': 'HMM file not found on the server. Please submit again.'})
                    if task_return[0] == 103:
                        return render(request, 'error.html', {'error_message': 'HMMSEARCH_RUNTIME_ERROR'})
                    context = {
                        'output_file_name': task_return[0],
                    }
                    return render(request, 'gene_family/gene_fam_identification_result.html', context)

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE_EXCEPTION'})
        # Render default view
        return render(request, 'gene_family/gene_hmmsearch_analysis.html')

    elif request.method == "POST":
        # Task submission
        # Obtain query from frontend form submission
        species_query = request.POST.get("query_species")  # by HTML tag name
        hmm_query = request.POST.get("hmm_query")
        sequence_evalue_threshold = request.POST.get("sequence_evalue")
        domain_evalue_threshold = request.POST.get("domain_evalue")
        try:
            # Validate if variable sequence_evalue_threshold and domain_evalue_threshold are numbers
            int(float(sequence_evalue_threshold)) / 2
            int(float(domain_evalue_threshold)) / 2
            # Make sure both numbers are positive
            if int(float(sequence_evalue_threshold)) < 1 or int(float(domain_evalue_threshold)) < 1:
                raise ValueError
        except ValueError:
            return render(request, 'error.html', {'error_message': 'INVALID_QUERY'})

        # Obtain filename for accessing background annotation
        species_filename = MainConfiguration.query_tables(str(species_query) + "_filename")

        # Deploy task
        # params: species_filename, hmm_query_filename, seq_evalue_threshold, dom_evalue_threshold
        async_result = celery_tasks.genefam_identification_hmmsearch_handler.delay(
            species_filename, hmm_query, int(float(sequence_evalue_threshold)), int(float(domain_evalue_threshold))
        )
        # Return async task manager, awaiting frontend recognition
        return redirect("/HEMUdb/asynctasks/onhold?taskid=%s&dir=%s" %
                        (async_result.task_id, "/HEMUdb/genefam/hmm"))


def genefam_expheatmapgen_async(request):
    if request.method == "GET":
        return render(request, 'tf_mainpage.html',
                      {
                          "zea_tissue_list": str(", ".join(TFHeatmapSampleConfiguration.return_sample("zea").keys())),
                      })

    elif request.method == "POST":  # Generate Heatmaps
        try:
            final_tf_fam_id = []
            final_tf_indv_id = []
            final_tissueid = []
            final_sampleid = []

            tf_query_list = request.POST.get('tf_query').split(';')
            sample_query_list = request.POST.get('sample_query').split(';')
            species_query = request.POST.get('species_query')

            # print(tf_query_list, sample_query_list, species_query)
            # Process TF query
            for indv_tf in tf_query_list:
                if str(indv_tf).startswith("fam:"):
                    final_tf_fam_id.append(str(indv_tf).lstrip("fam:"))
                elif str(indv_tf).startswith("id:"):
                    final_tf_indv_id.append(str(indv_tf).lstrip("id:"))
            final_TF_fam_list = TFHeatmapDataCurator.TF_fam_to_geneid_query(species_query, final_tf_fam_id,
                                                                            final_tf_indv_id)
            # final_TF_fam_list: [[TF_geneid, TF_family], [TF_geneid, TF_family], ...]

            # Process sample query
            for indv_sample in sample_query_list:
                if str(indv_sample).startswith("tis:"):
                    final_tissueid.append(str(indv_sample).lstrip("tis:"))
                elif str(indv_sample).startswith("sp:"):
                    final_sampleid.append(str(indv_sample).lstrip("sp:"))

            final_heatmap_df = TFHeatmapDataCurator.TF_heatmap_df_builder(species_query, final_TF_fam_list,
                                                                          final_tissueid)
            heatmap_filename = TF_heatmap_generator.TF_tpm_heatmap(final_heatmap_df)

            return render(request, 'tf_heatmap_display.html',
                          {
                              'heatmap_filename': heatmap_filename,
                              'species_query': species_query,
                              # 'confirmation_msg': confirmation_msg,
                          })
        except:
            return HttpResponse("Some error occurred, please try again.")


# ==============MODULE: TE==============
# 4.1 TE HOMEPAGE RENDERER
def te_exp_init(request):
    global exp_df

    if request.method == "GET":
        return render(request, 'TE/TE_search_main.html')
    elif request.method == "POST":
        species_by_accession = request.POST.get('species_query_accession')
        species_by_teid = request.POST.get('species_query_teid')

        if species_by_accession:  # by-sample-accession search
            te_query_accession_list = request.POST.get('te_query_accession').split(";")
            for te_query_accession in te_query_accession_list:
                if species_by_accession == "coix":
                    exp_df = TEBySampleDataCurator.TE_exp_df_builder_sample(te_query_accession,
                                                                            MainConfiguration.query_tables('coix_te'),
                                                                            # Expression matrix
                                                                            )  # Detailed sample info
                elif species_by_accession == "zea":
                    exp_df = TEBySampleDataCurator.TE_exp_df_builder_sample(te_query_accession,
                                                                            MainConfiguration.query_tables('zea_te'),
                                                                            # Expression matrix
                                                                            )  # Detailed sample info
                elif species_by_accession == "sorghum":
                    exp_df = TEBySampleDataCurator.TE_exp_df_builder_sample(te_query_accession,
                                                                            MainConfiguration.query_tables(
                                                                                'sorghum_te'),  # Expression matrix
                                                                            )  # Detailed sample info
                else:
                    _error_message = "illegal query"
                TE_expression_plt_generator.TE_bysample_plt(exp_df, te_query_accession)

            return render(request, 'TE/TE_bysample_result.html',
                          {
                              'query_list': te_query_accession_list
                          })

        elif species_by_teid:  # by-TE family ID search
            te_query_TEid_list = request.POST.get('te_query_teid').split(";")
            for te_query_TEid in te_query_TEid_list:
                if species_by_teid == "coix":
                    TE_class, TE_class_group, exp_df = TEByFamilyDataCurator.TE_exp_df_builder_family(te_query_TEid,
                                                                                                      MainConfiguration.query_tables(
                                                                                                          'coix_te'),
                                                                                                      # Expression matrix
                                                                                                      )  # Detailed sample info
                elif species_by_teid == "zea":
                    TE_class, TE_class_group, exp_df = TEByFamilyDataCurator.TE_exp_df_builder_family(te_query_TEid,
                                                                                                      MainConfiguration.query_tables(
                                                                                                          'zea_te'),
                                                                                                      # Expression matrix
                                                                                                      )  # Detailed sample info
                elif species_by_teid == "sorghum":
                    TE_class, TE_class_group, exp_df = TEByFamilyDataCurator.TE_exp_df_builder_family(te_query_TEid,
                                                                                                      MainConfiguration.query_tables(
                                                                                                          'sorghum_te'),
                                                                                                      # Expression matrix
                                                                                                      )  # Detailed sample info
                else:
                    _error_message = "illegal query"
                TE_expression_plt_generator.TE_byfamily_plt(exp_df, te_query_TEid)

            return render(request, 'TE/TE_byfamily_result.html',
                          {
                              'query_list': te_query_TEid_list
                          })

        else:
            return HttpResponse("Legal query not detected.")


def te_insertion(request):
    if request.method == "GET":
        return render(request, 'TE/TE_insertion.html')


# ==============MODULE: DATA BROWSER - JBrowse2==============
# 5.1 JBROWSE2 RENDERER
def jbrowse_catalog_render(request):
    if request.method == "GET":
        # Accession ID, jbrowse/?id=GCF_000005005.1
        full_name = jbrowse_species_catalog.jbrowse_config_query(str(request.GET.get('id')))
        if full_name:  # Determine whether client is sending a query
            if not full_name == -1:  # Capture key error
                return render(request, 'Jbrowse/jbrowse_insession.html',  # Sending a valid query
                              {'species_full_name': full_name})
        return render(request, 'Jbrowse/jbrowse_main_catalog.html')  # Sending an invalid query or not sending a query


# ==============MODULE: BLAST - SequenceServer==============
def seqserver(request):
    if request.method == "GET":
        return render(request, 'BLAST/sequenceserver_display.html', {'seq_transfer': ""})  # No autofill


# ==============MODULE: SUPP==============
# 6.1 USER GUIDE RENDERER
def user_guide(request):
    return render(request, 'user_guide.html')

# ==============MODULE: MULTI-OMICS DATA ANALYSIS==============
# 2.1 miRNA-seq
def miRNA_expression_query_async(request):
    if request.method == 'GET':
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Task completed, return [query_list, query_list_full, query_format, query_species]
                    context = {
                        'query_list': task_return[0],
                        'query_list_full': task_return[1],
                        'query_format': task_return[2],
                        'query_species': task_return[3],
                        'last_query': ";".join(task_return[0]),
                    }
                    return render(request, 'miRNA/miRNA_search_main.html', context)

                # task_return == None means no database queue hit, probably due to faulty input.
                return render(request, 'miRNA/miRNA_search_main.html',
                              {'error_message': 'Invalid entry. Please check gene nomenclature and submit the query '
                                                'again.'})

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE'})
        # Render default view
        return render(request, 'miRNA/miRNA_search_main.html')

    elif request.method == "POST":
        # Task submission
        # Obtain query from frontend form submission
        query_list = request.POST.get('main_query').split(';')  # Original query data
        query_species = request.POST.get('query_species')  # Species to query
        query_format = request.POST.get('optionsRadios')  # Counts/Nomalized

        # Obtain sql tables
        query_exp_table = MainConfiguration.query_tables(str(query_species) + "_miRNA_exp")

        # Deploy task
        async_result = celery_tasks.gene_expression_query_handler.delay(
            query_species,
            query_exp_table,
            query_list,
            query_format,
        )
        # Return page containing task id, awaiting frontend recognition
        return render(request, 'miRNA/miRNA_search_main.html', {'task_id': async_result.task_id})

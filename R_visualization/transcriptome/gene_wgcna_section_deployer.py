#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project ：HEMU_Database_Main 
@File    ：gene_wgcna_section_deployer.py
@Author  ：Edward Zhu, zhuyzh37@mail2.sysu.edu.cn
@Date    ：2023/03/04 14:43 
@IDE     ：PyCharm 
-----------------------------
Description:
Input:
Output:

One-line usage note: python3 gene_wgcna_section_deployer.py ,
- ARG1:

Environment requirement (if any): base
- Dependencies: None
-----------------------------
Notes:

-----------------------------
File Revisions:
    2023/03/04: Version 1 - Creation
"""
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import random
import os

# Used for passing pandas dataframe variable to R
pandas2ri.activate()  # Used for variable conversion (pd.DataFrame)


def wgcna_section1_deployer(sampleGeneData, RcCutoff, samplePerc, GeneNum, cutmethod, rscut, datatype, anamethod,
                            dirname):

    # --- Read R script ---
    rscript_fh = open("Mainapp/R_scripts/wgcna_plotter.R")
    # readlines() generates a list, concat pending.
    rscript = "".join(rscript_fh.readlines())
    rscript_fh.close()
    robjects.r(rscript)
    # --- Finish reading R script ---

    robjects.globalenv['sampleGeneData'] = sampleGeneData
    robjects.globalenv['RcCutoff'] = RcCutoff
    robjects.globalenv['samplePerc'] = samplePerc
    robjects.globalenv['GeneNum'] = GeneNum
    robjects.globalenv['cutmethod'] = cutmethod
    robjects.globalenv['rscut'] = rscut
    robjects.globalenv['datatype'] = datatype
    robjects.globalenv['anamethod'] = anamethod
    robjects.globalenv['dirname'] = dirname

    # Execute function
    robjects.r['gene_wgcna_step1'](sampleGeneData, RcCutoff, samplePerc, GeneNum, cutmethod, rscut, datatype, anamethod,
                                   dirname)

    return dirname


def wgcna_section2_deployer(sftPower, minModuleSize, mergeCutHeight, dirname):

    # Detect and delete existing plots
    if not os.path.exists('Mainapp/static/Temp_R_wgcna/' + dirname):
        return RuntimeError("Project not found.")

    # --- Read R script ---
    rscript_fh = open("Mainapp/R_scripts/wgcna_plotter.R")
    # readlines() generates a list, concat pending.
    rscript = "".join(rscript_fh.readlines())
    rscript_fh.close()
    robjects.r(rscript)
    # --- Finish reading R script ---

    robjects.globalenv['sftPower'] = sftPower
    robjects.globalenv['minModuleSize'] = minModuleSize
    robjects.globalenv['mergeCutHeight'] = mergeCutHeight
    robjects.globalenv['dirname'] = dirname

    # Execute function
    robjects.r['gene_wgcna_step2'](sftPower, minModuleSize, mergeCutHeight,
                                   dirname)

    return dirname


def wgcna_section3_deployer(traitData, dirname):

    # Detect and delete existing plots
    if not os.path.exists('Mainapp/static/Temp_R_wgcna/' + dirname):
        return RuntimeError("Project not found.")

    # --- Read R script ---
    rscript_fh = open("Mainapp/R_scripts/wgcna_plotter.R")
    # readlines() generates a list, concat pending.
    rscript = "".join(rscript_fh.readlines())
    rscript_fh.close()
    robjects.r(rscript)
    # --- Finish reading R script ---

    robjects.globalenv['traitData'] = "Mainapp/static/Temp_R_genefam/" + dirname + "/" + traitData
    robjects.globalenv['dirname'] = dirname

    # Execute function
    robjects.r['gene_wgcna_step3'](traitData, dirname)

    return dirname

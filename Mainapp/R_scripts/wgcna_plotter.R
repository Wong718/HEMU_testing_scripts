
gene_wgcna_step1 <- function (sampleGeneData,
                              RcCutoff, samplePerc, GeneNum, cutmethod, rscut, datatype, anamethod, dirname){
  # Data cleaning and sft calculation
  # ====================================
  # 'datatype''Data type: count or FPKM. count means readcount, FPKM represent normalized count,eg: FPKM TPM CPM RPKM,
  # 'anamethod''Data transformat method: For count: varianceStabilizingTransformation or cpm; For FPKM, rawFPKM or logFPKM,
  # 'RcCutoff''Noise cutoff: count/normalized count lower than ... was thought to be noise. Default: count=>10, Normalized count=>1',
  # 'samplePerc''Sample percentage: parameter for noise remove. xx percent of all samples have readcount/FPKM > cutoff. Default: 0.3',
  # 'GeneNum''Gene number for WGCNA: how many gene you want retained for WGCNA after noise remove.',
  # 'cutmethod''MAD or SVR, Default: MAD',
  # 'rscut''sft Power cutoff: Power cutoff',
  
  #setwd(paste0("global-temp/", dirname))
  rm(list=ls())
  
  # Check if dependencies are fully installed
  # suppressMessages(suppressWarnings(if (!require('logr')) install.packages('logr')));
  suppressMessages(suppressWarnings(if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")))
  suppressMessages(suppressWarnings(if (!require('devtools')) install.packages('devtools')));
  suppressMessages(suppressWarnings(if (!require('WGCNA')) BiocManager::install('WGCNA',update = FALSE)));
  suppressMessages(suppressWarnings(if (!require('tidyverse')) install.packages('tidyverse')));
  suppressMessages(suppressWarnings(if (!require('ShinyWGCNA')) devtools::install_github("ShawnWx2019/WGCNAShinyFun",ref = "master")));
  suppressMessages(suppressWarnings(if (!require('DESeq2')) BiocManager::install('DESeq2',update = FALSE)));
  suppressMessages(suppressWarnings(if (!require('ggprism')) BiocManager::install('ggprism',update = FALSE)));
  suppressMessages(suppressWarnings(if (!require('patchwork')) BiocManager::install('patchwork',update = FALSE)));
  # Enable multi-threading
  enableWGCNAThreads(nThreads = 2)
  # For macOS
  #allowWGCNAThreads()
  
  # Testing parameters
  # sampleGeneData = "WGCNA_fpkm_data.csv"
  # RcCutoff = 1
  # samplePerc = 0.9
  # GeneNum = 2000
  # cutmethod = "MAD"
  # rscut = 0.8
  # datatype = "FPKM"
  # anamethod = "rawFPKM"

  # Conversion from character (python-relay) to number
  RcCutoff = as.numeric(RcCutoff)
  samplePerc = as.numeric(samplePerc)
  GeneNum = as.numeric(GeneNum)
  rscut = as.numeric(rscut)
  
  dx = read.table(sampleGeneData, header=T, sep = ",") # rowname-genes, colname-samples
  # print(summary(dx))
  ngenes = nrow(dx) # Obtain total number of genes
  GeneNumCut = 1-(GeneNum/ngenes) # Obtain proportion of cut threshold
  
  # Dataframe format conversion, as.character -> as.numeric
  vec_datexpr_rownames = dx$X
  dx = dx[,-1]
  dx = as.data.frame(lapply(dx,as.numeric))
  dx = cbind(vec_datexpr_rownames, dx)
  rownames(dx) = vec_datexpr_rownames

  # --NOISE REMOVAL--
  datExpr <- getdatExpr(rawdata = dx, RcCutoff = RcCutoff, samplePerc = samplePerc, datatype = datatype, method = anamethod)
  datExpr = datExpr[,-1]

  # --GENE FILTRATION--
  datExpr <- getdatExpr2(datExpr = datExpr, GeneNumCut = GeneNumCut, cutmethod = cutmethod)

  # --SFT PREDICTION--
  step1_sft = getpower(datExpr = datExpr, rscut = rscut)
  # Save scale-independence and mean connectivty plot
  png(filename = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step1-sft.png"), width=800, height=400)
  step1_sft$plot
  dev.off()
  write_csv(step1_sft$sft, file = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step1-sft-data.csv"))
  print(paste0( "SFT prediction completed, The recommanded power is: ", step1_sft$power))

  # Save workspace
  save(datExpr, RcCutoff, samplePerc, GeneNum, cutmethod, rscut, datatype, anamethod,
       file=paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step1-vars.RData"))
}

gene_wgcna_step2 <- function (sftPower, minModuleSize, mergeCutHeight, dirname){
  # Gene expression network construction
  # ====================================
  # 'sftPower''SFT power.',
  # 'minModuleSize''Minimal module size: The gene number of minimal module. Default: 30',
  # 'mergeCutHeight''Merge cuttree height: tree height lower than this value will be merged. Default: 0.25',
  
  #setwd(paste0("global-temp/", dirname))
  rm(list=ls())
  
  # Check if dependencies are fully installed
  suppressMessages(suppressWarnings(library(WGCNA)));
  suppressMessages(suppressWarnings(library(ShinyWGCNA)));
  suppressMessages(suppressWarnings(library(tidyverse)));
  # Enable multi-threading
  enableWGCNAThreads(nThreads = 2)
  
  # Testing parameters
  # sftPower = 6
  # minModuleSize = 30
  # mergeCutHeight = 0.25
  
  # Load data from previous steps
  load(paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step1-vars.RData"))
  # Prepare log destination, open file handler
  tmp <- file.path(paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step2-network-construction.log"))
  lf <- log_open(tmp)
  
  # Set max block size
  maxBlocksize = ncol(datExpr)
  # 'maxBlocksize''max block size: For block-wised network construction method, the block size, which is set to be all filtered genes.'
  
  log_print("--GENE EXPRESSION NETWORK CONSTRUCTION STARTED--")
  tryCatch({
    step2_network = getnetwork(datExpr = datExpr, power = sftPower, 
                               minModuleSize = minModuleSize, mergeCutHeight = mergeCutHeight, maxBlocksize = maxBlocksize)
    # Visualization using the hierarchical clustering dendrogram
    gene_tree = step2_network$net$dendrograms[[1]]
    png(filename = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step2-network-dendro.png"), width=1000, height=600, res=100)
    plotDendroAndColors(
      gene_tree,
      step2_network$moduleColors[step2_network$net$blockGenes[[1]]],
      dendroLabels = F
    )
    dev.off()
    write_csv(step2_network$MEs_col, file = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step2-module-eigengenes.csv"))
    write_csv(step2_network$Gene2module, file = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step2-module-allgenes.csv"))
    
    log_print("Network construction completed.")
  },error = function(e) {
    log_print("Network construction terminated with non-zero exit status, inspect log for more details.")
    q(status = 1)
  })
  
  # Export parameters to log
  params = data.frame(
    param = c("datatype","anamethod","RcCutoff","samplePerc","GeneNum","cutmethod","rscut","sftPower","minModuleSize","mergeCutHeight","maxBlocksize"),
    value = c(datatype,anamethod,RcCutoff,samplePerc,GeneNum,cutmethod,rscut,sftPower,minModuleSize,mergeCutHeight,maxBlocksize)
  )
  log_print("Parameters used in this step: \n")
  log_print(params)
  
  # Save workspace
  #save.image("step2.RData")
  save(step2_network, datExpr, 
       RcCutoff, samplePerc, GeneNum, cutmethod, rscut, datatype, anamethod,
       sftPower, minModuleSize, mergeCutHeight, maxBlocksize,
       file=paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step2-vars.RData"))
  
  # Close log file handle
  log_close()
  writeLines(readLines(lf))
}

gene_wgcna_step3 <- function (traitData, dirname){
  # Module-trait correlation analysis
  # ====================================
  # 'traitData''trait-sample information, trait represented with 0 or 1',
  
  #setwd(paste0("global-temp/", dirname))

  suppressMessages(suppressWarnings(library(ggplot2)));
  suppressMessages(suppressWarnings(library(WGCNA)));
  suppressMessages(suppressWarnings(library(ShinyWGCNA)));
  suppressMessages(suppressWarnings(library(tidyverse)));
  suppressMessages(suppressWarnings(library(ggprism)));
  suppressMessages(suppressWarnings(library(patchwork)));
  suppressMessages(suppressWarnings(if (!require('ComplexHeatmap')) BiocManager::install('ComplexHeatmap',update = FALSE)));
  # Enable multi-threading
  enableWGCNAThreads(nThreads = 2)
  
  # Testing parameters
  # traitData = "WGCNA_data_trait.csv"
  
  # Load data from previous steps
  load(paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step2-vars.RData"))
  # set options for scientific number notation
  options(scipen = 6)
  # Prepare log destination, open file handler
  tmp <- file.path(paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step3-module-trait-analysis.log"))
  lf <- log_open(tmp)
  
  log_print("--MODULE-TRAIT CORRELATION ANALYSIS STARTED--")
  phen <- read.csv(traitData)
  
  # Bind sample-id column to the dataframe index
  colnames(phen)[1] = "sample_id"
  phen <- left_join(
    data.frame(
      sample_id = rownames(datExpr)
    ),phen,"sample_id"
  ) %>% 
    column_to_rownames("sample_id")
  
  tryCatch({
    step3_mdule_trait <- getMt(
      phenotype = phen,
      nSamples = nrow(datExpr),
      moduleColors = step2_network$moduleColors,
      datExpr = datExpr
    )
    log_print("Module-trait correlation analysis completed.")
  },error = function(e) {
    log_print("Module-trait correlation analysis terminated with non-zero exit status, inspect log for more details.")
    q(status = 1)
  })
  
  mod_color = gsub(pattern = "^..",replacement = "",rownames(step3_mdule_trait$modTraitCor))
  mod_color_anno = setNames(mod_color,rownames(step3_mdule_trait$modTraitCor))
  
  left_anno = rowAnnotation(
    Module = rownames(step3_mdule_trait$modTraitCor),
    col = list(Module = mod_color_anno),
    show_legend = F,
    show_annotation_name = F
  )
  
  # Generate data matrix, generate heatmap
  log_print("Generating heatmap")
  htmap_df = as.matrix(step3_mdule_trait$modTraitCor)
  
  # Save heatmap to static image
  png(filename = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step3-module-trait-heatmap.png"), width=2000, height=2000, res=200)
  pheatmap(htmap_df,
           display_numbers = TRUE,
           color = colorRampPalette(c("#00AFBB", "white", "#c93756"))(256),
           fontsize=10, border_color = "grey60",
           main = "Module-trait heatmap",
           treeheight_row=50, treeheight_col = 30,
           cellwidth = 50, cellheight = 40
  )
  dev.off()
  write_csv(as.data.frame(htmap_df), file = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step3-module-trait-correlation.csv"))
  
  # Export parameters to log
  params = data.frame(
    param = c("datatype","anamethod","RcCutoff","samplePerc","GeneNum","cutmethod","rscut","sftPower","minModuleSize","mergeCutHeight","maxBlocksize"),
    value = c(datatype,anamethod,RcCutoff,samplePerc,GeneNum,cutmethod,rscut,sftPower,minModuleSize,mergeCutHeight,maxBlocksize)
  )
  log_print("Parameters used in this step: \n")
  log_print(params)
  
  log_print("--ALL WORK FINISHED, CHECK RESULTS--")
}

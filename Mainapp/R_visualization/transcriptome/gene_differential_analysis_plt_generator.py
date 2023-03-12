import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import os
import random

# Used for passing pandas dataframe variable to R
pandas2ri.activate()  # Used for variable conversion (pd.DataFrame)


def GeneDifferentialAnalysis(DE_data_raw, DE_group_list, DE_group_color_list,
                             logfc_threshold, pvalue_threshold, heatmap_gene_count,
                             group1_name, group2_name):
    # Generate 7-digit random number for folder nomination
    output_folder_name = str(random.randint(1000000000, 9999999999))

    os.makedirs('Mainapp/static/Temp_R_DEprojects/' + output_folder_name, exist_ok=True)
    # print("Folder made done, %s" % output_folder_name)

    robjects.globalenv['DE_data_raw'] = DE_data_raw
    robjects.globalenv['DE_group_list'] = DE_group_list
    robjects.globalenv['DE_group_color_list'] = DE_group_color_list
    robjects.globalenv['logfc_threshold'] = logfc_threshold
    robjects.globalenv['pvalue_threshold'] = pvalue_threshold
    robjects.globalenv['heatmap_gene_count'] = heatmap_gene_count
    robjects.globalenv['group1_name'] = group1_name
    robjects.globalenv['group2_name'] = group2_name
    robjects.globalenv['output_folder_name'] = output_folder_name

    rscript_DE_analysis = '''

    DE_analysis <- function (){
      library(limma)
      library(tidyverse)
      library(plotly)
      library(reshape2)
      
      library(FactoMineR) # PCA
      library(factoextra)
      
      library(pheatmap) # Heatmap
      library(ggcorrplot) # Correlation plot
      
      library(viridis)
      
      #print("Library Load Complete.")
      heatmap_top_gene_count = heatmap_gene_count
      
      #print(summary(DE_data_raw))
      #print(rownames(DE_data_raw))
      #rownames(DE_data_raw) = DE_data_raw[,1]
      #DE_data_raw = DE_data_raw[,-1]
      
      DE_group_list = as.character(DE_group_list)
      group_list = factor(DE_group_list)
      group_color_list = as.character(DE_group_color_list)
      
      exprSet = DE_data_raw
    
      # Normalization of TPM data
      exprSet_boxplt_data = melt(data=exprSet)
      colnames(exprSet_boxplt_data) = c("sample_id", "Raw_TPM")
      
      plt_boxplot_beforenorm = ggplot(data=exprSet_boxplt_data) + 
        geom_boxplot(aes(x=sample_id, y=Raw_TPM), outlier.shape = NA, alpha=0.7,
                     outlier.size=0.1,
                     fill = group_color_list) +
        coord_cartesian(ylim = c(0, 60)) +
        theme(panel.grid = element_blank(), 
              panel.background = element_blank(), 
              axis.line = element_line(colour = 'black'),
              axis.text.x = element_text(size = 10, hjust = 1, angle = 45),
              axis.text.y = element_text(size = 10, hjust = 1),
              axis.title.x = element_blank()
        )
      
      # Normalize Between Arrays
      exprSet=normalizeBetweenArrays(exprSet)
      
      exprSet_boxplt_data = melt(data=exprSet)[,-1]
      colnames(exprSet_boxplt_data) = c("sample_id", "Normalized_TPM")
      
      plt_boxplot_afternorm = ggplot(data=exprSet_boxplt_data) + 
        geom_boxplot(aes(x=sample_id, y=Normalized_TPM), outlier.shape = NA, alpha=0.7,
                     outlier.size=0.1,
                     fill = group_color_list) +
        coord_cartesian(ylim = c(0, 60)) +
        theme(panel.grid = element_blank(), 
              panel.background = element_blank(), 
              axis.line = element_line(colour = 'black'),
              axis.text.x = element_text(size = 10, hjust = 1, angle = 45),
              axis.text.y = element_text(size = 10, hjust = 1),
              axis.title.x = element_blank()
        )
      
      #print("Normalization of TPM data complete.")
      exprSet <- log2(exprSet+1)
      
      # Execution of DE analysis
      dat <- exprSet
      design=model.matrix(~factor( group_list ))
      fit=lmFit(dat,design)
      fit=eBayes(fit)
      options(digits = 4)
      topTable(fit,coef=2,adjust='BH')
      bp=function(g){
        
        df=data.frame(gene=g,stage=group_list)
        p <- ggboxplot(df, x = "stage", y = "gene",
                       color = "stage", palette = "jco",
                       add = "jitter")
        #  Add p-value
        p + stat_compare_means()
      }
      deg=topTable(fit,coef=2,adjust='BH',number = Inf)
      
      #print("Execution completed of DE, step 1")
      if(T){
        logFC_t=2
        deg$status=ifelse(deg$P.Value>0.05,'NS',
                          ifelse( deg$logFC > logFC_t,'UP',
                                  ifelse( deg$logFC < -logFC_t,'DOWN','NS') )
        )
        deg$symbol=rownames(deg)
      }
      
      ### Output DEG results into a single .csv file.
      write.csv(deg, paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/differential_gene.csv"), row.names = T)
      
      ### Density plot of expression ###
      DEraw_density_data = melt(data=DE_data_raw)
      colnames(DEraw_density_data) = c("sample_id", "Raw_TPM")
      plt_sample_exp_density = ggplot(data=DEraw_density_data) +
          geom_density(aes(x=Raw_TPM, fill=sample_id), adjust=1.5) +
          coord_cartesian(xlim = c(0, 300)) +
          scale_fill_viridis(option = 'cividis', discrete = T) +
          theme_minimal() +
          facet_wrap(~sample_id) +
          theme(#panel.grid = element_blank(), 
                #panel.background = element_blank(), 
                axis.line = element_line(colour = 'black'),
                axis.text.x = element_text(size = 10, hjust = 1, angle = 45),
                axis.text.y = element_text(size = 10, hjust = 1),
                legend.position = 'none',
          ) +
          labs(y="Density", x="TPM (Before Normalization)")
      
            
      ### Sample Correlation Plot ###
      # Use TPM values after both limma normalization and log2(N+1) transformation
      exp_correlation = cor(exprSet)
      exp_correlation_p <- cor_pmat(exprSet)
    
      plt_sample_correlation = ggcorrplot(
        exp_correlation,
        outline.color = "black",
        #method = "circle",
        ggtheme = ggplot2::theme_bw,
        lab_size = 4, p.mat = exp_correlation_p, sig.level = 0.05, insig = "blank", pch = 4, pch.col = "black", pch.cex = 5, tl.cex = 12,
        colors = c("#839EDB", "white", "#FF8D8D"),lab = T, title=paste0("Correlation: ", group1_name, " vs ", group2_name)) +
      theme(
        axis.text.x = element_text(size = 10, hjust = 1),
        axis.text.y = element_text(size = 10, hjust = 1),
      )
      
      #print("Execution completed of DE, step 2. Generating volcano plot")
      plt_volcano = ggplot(data=deg, aes(x=logFC, y=-log10(P.Value), color=status)) +
        geom_point(size=0.6, alpha=0.6) +
        scale_color_manual(values =c("DOWN"='#00AFBB',"NS"='#eaeaea',"UP"="#f38181")) +
        
        geom_vline(aes(xintercept=logFC_t),color="#40514e",size=0.5,linetype="dashed") +
        geom_vline(aes(xintercept=-logFC_t),color="#40514e",size=0.5,linetype="dashed") +
        geom_hline(aes(yintercept=-log10(as.numeric(pvalue_threshold))),color="#40514e",size=0.5,linetype="dashed") + # -log10(0.05)
        
        theme(panel.grid = element_blank(), 
              panel.background = element_blank(), 
              axis.line = element_line(colour = 'black'),
              axis.text.x = element_text(size = 10, hjust = 1),
              axis.text.y = element_text(size = 10, hjust = 1),
              #axis.title.x = element_blank(),
              #legend.position = 'none',
        ) +
        labs(title=paste0("Volcano: ", group1_name, " vs ", group2_name))
      
      ### PCA Analysis ###
      #print("Performing PCA analysis")
      raw_PCA_data = t(exprSet)
      dat.pca = PCA(raw_PCA_data, graph = F)
      
      plt_pca_samples = fviz_pca_ind(dat.pca,
                                     #geom.ind = "point", # show points only (but not "text")
                                     col.ind = DE_group_list, # color by groups
                                     labelsize = 3,
                                     palette = c("#00AFBB", "#E7B800"),
                                     addEllipses = TRUE, # Concentration ellipses 
                                     legend.title = "Groups") +
        
        theme(panel.grid = element_blank(), 
              panel.background = element_blank(), 
              axis.line = element_line(colour = 'black'),
              axis.text.x = element_text(size = 10, hjust = 1),
              axis.text.y = element_text(size = 10, hjust = 1),
              #legend.position = 'none',
        ) +
        labs(title=paste0("PCA-samples: ", group1_name, " vs ", group2_name))
      
      ### Heatmap of top most DE genes ###
      #print("Plotting heatmap")
      
      DE_data_heatmap = deg[which(deg$status %in% c("UP", "DOWN")),]
      DE_data_heatmap$logFC_abs = abs(DE_data_heatmap$logFC)
      DE_data_heatmap = DE_data_heatmap[order(-DE_data_heatmap$logFC_abs),]
      
      DE_data_heatmap_toplist = DE_data_heatmap[1:heatmap_top_gene_count,8]
      
      TPM_data_raw_heatmap = DE_data_raw
      TPM_data_raw_heatmap$symbol = rownames(TPM_data_raw_heatmap)
      top_DE_rownames = rownames(TPM_data_raw_heatmap)
      
      TPM_data_top_heatmap = TPM_data_raw_heatmap[which(TPM_data_raw_heatmap$symbol %in% DE_data_heatmap_toplist),-ncol(TPM_data_raw_heatmap)]
      TPM_data_top_heatmap = as.matrix(TPM_data_top_heatmap)
      #TPM_data_top_heatmap = apply(TPM_data_top_heatmap,2,as.numeric)
      
      TPM_data_top_heatmap[TPM_data_top_heatmap==0] = NA
      TPM_data_top_heatmap[is.na(TPM_data_top_heatmap)] = min(TPM_data_top_heatmap,na.rm = T)*0.01
      
      #TPM_data_top_heatmap[is.na(TPM_data_top_heatmap)] = 0.01
      TPM_data_top_heatmap = TPM_data_top_heatmap[apply(TPM_data_top_heatmap, 1, function(x) sd(x)!=0),]
      # Eliminate pheatmap error: NA/NaN/Inf in foreign function call
      
      rownames(TPM_data_raw_heatmap) = top_DE_rownames
      
      plt_heatmap_topDE = pheatmap(log10(TPM_data_top_heatmap),
                                   display_numbers = TRUE,
                                   #color = colorRampPalette(c("#00AFBB", "white", "#c93756"))(256),
                                   color = viridis(256, option='cividis'),
                                   fontsize=7, border_color = "grey60",
                                   main = paste0("Heatmap - top DE genes: ", group1_name, " vs ", group2_name, "\n(log10 transformed)"),
                                   treeheight_row=50, treeheight_col = 30,
                                   cellwidth = 25, cellheight = 8,
                                   silent=T,
      )
      
      ggsave(paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/boxplot_beforenorm.png"), 
             plot=plt_boxplot_beforenorm, dpi = 400, width=5, height=5)
      output_box_beforenorm <- plotly::ggplotly(plt_boxplot_beforenorm)
      htmlwidgets::saveWidget(output_box_beforenorm, paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/boxplot_beforenorm.html"),
        selfcontained = TRUE)
        
      ggsave(paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/boxplot_afternorm.png"), 
           plot=plt_boxplot_afternorm, dpi = 400, width=5, height=5)
      output_box_afternorm <- plotly::ggplotly(plt_boxplot_afternorm)
      htmlwidgets::saveWidget(output_box_afternorm, paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/boxplot_afternorm.html"),
        selfcontained = TRUE)
        
      
      ggsave(paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/exp_density.png"), 
           plot=plt_sample_exp_density, dpi = 400, width=5, height=5)
      output_sample_exp_density <- plotly::ggplotly(plt_sample_exp_density)
      htmlwidgets::saveWidget(output_sample_exp_density, paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/exp_density.html"),
        selfcontained = TRUE)
        
    
      ggsave(paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/exp_correlation.png"), 
           plot=plt_sample_correlation, dpi = 400, width=5, height=5)
      output_sample_correlation <- plotly::ggplotly(plt_sample_correlation)
      htmlwidgets::saveWidget(output_sample_correlation, paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/exp_correlation.html"),
        selfcontained = TRUE)
      
      ggsave(paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/volcano.png"), 
        plot=plt_volcano, dpi = 400, width=7, height=5)
      output_volcano <- plotly::ggplotly(plt_volcano)
      htmlwidgets::saveWidget(output_volcano, paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/volcano.html"),
        selfcontained = TRUE)
      

      output_pca_samples <- plotly::ggplotly(plt_pca_samples)
      htmlwidgets::saveWidget(output_pca_samples, paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/pca_samples.html"),
        selfcontained = TRUE)
      ggsave(paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/pca_samples.png"), 
        plot=plt_pca_samples, dpi = 400, width=7, height=5)
      
      ggsave(paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/topDE_heatmap.png"), 
             plot=plt_heatmap_topDE, dpi = 400,
             width = length(DE_group_list), height = as.numeric(heatmap_top_gene_count)/7)

    }
    '''

    robjects.r(rscript_DE_analysis)
    robjects.r['DE_analysis']()  # Execute function

    return output_folder_name

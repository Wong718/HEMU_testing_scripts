import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

pandas2ri.activate()  # Used for variable conversion (pd.DataFrame)


def TE_bysample_plt(_te_df, sample_id):

    robjects.globalenv['te_dataframe'] = _te_df
    robjects.globalenv['filename_id'] = sample_id

    plt_rscript = '''

    TE_sample_plotter <- function (te_dataframe, filename_id){
        rm(list=ls())

        library(tidyverse)
        #library(patchwork)
        library(plotly)
        te_total_expressed = te_dataframe
        output_id = filename_id
        
        te_total_expressed = te_total_expressed[,-1]
        te_total_expressed$TE_class = as.factor(te_total_expressed$TE_class)
        te_total_expressed$TE_class_group = as.factor(te_total_expressed$TE_class_group)
        te_total_expressed$fpkm = as.numeric(te_total_expressed$fpkm)
        te_total_expressed$tpm = as.numeric(te_total_expressed$tpm)
        
        
        plt_TE_supfam = ggplot(data=te_total_expressed,
                  aes(x=TE_class, y=tpm)) +
          geom_jitter(aes(color=TE_class_group), width=0.2, size=0.5, alpha=0.4) +
          geom_boxplot(aes(fill=TE_class_group), alpha=0.6, notch = T, outlier.shape = NA) +
          
          theme(axis.text.x = element_text(angle = 0, size = 8),) +
          theme_bw() +
          scale_fill_brewer(palette = "BrBG") +
          scale_color_brewer(palette = "BrBG") +
          coord_flip() +
          theme(
                #panel.grid = element_blank(), 
                #panel.background = element_blank(), 
                axis.line = element_line(colour = 'black'),
                axis.text.x = element_text(size = 10, hjust = 1, angle = 45),
                axis.text.y = element_text(size = 10, hjust = 1),
            ) +
          ylim(1,200) +
          #facet_grid(TE_class_group~., drop = T) +
          xlab("TE Superfamily") + ylab("TPM") +
          labs(color = "Group", fill = "Group")
        
        
        plt_TE_group = ggplot(data=te_total_expressed, aes(x=TE_class_group, y=tpm)) +
          geom_jitter(aes(color=TE_class_group), width=0.1, size=0.5, alpha=0.2) +
          geom_boxplot(aes(fill=TE_class_group), alpha=0.6, notch = T, outlier.shape = NA) +
          scale_fill_brewer(palette = "BrBG") +
          scale_color_brewer(palette = "BrBG") +
          #coord_flip() +
          
          theme(
                panel.grid = element_blank(), 
                panel.background = element_blank(), 
                axis.line = element_line(colour = 'black'),
                axis.text.x = element_text(size = 10, hjust = 1, angle = 60),
                axis.text.y = element_text(size = 10, hjust = 1),
                legend.position = 'none'
            ) +
          ylim(1,200) +
          #facet_grid(TE_class_group~., drop = T) +
          xlab("TE Group") + ylab("TPM")
        
        
        output_plt_TE_supfam <- plotly::ggplotly(plt_TE_supfam)
        htmlwidgets::saveWidget(output_plt_TE_supfam, 
                                paste0("Mainapp/static/Temp_R_TE/", output_id,
                                       "_plt_TE_superfamily.html"),
                                selfcontained = TRUE)
        
        output_plt_TE_group <- plotly::ggplotly(plt_TE_group)
        htmlwidgets::saveWidget(output_plt_TE_group, 
                                paste0("Mainapp/static/Temp_R_TE/", output_id,
                                       "_plt_TE_group.html"),
                                selfcontained = TRUE)
    
    }
        
    '''

    robjects.r(plt_rscript)
    robjects.r['TE_sample_plotter'](_te_df, sample_id)  # Execute function


def TE_byfamily_plt(_te_df, TE_id):
    robjects.globalenv['te_dataframe'] = _te_df
    robjects.globalenv['filename_id'] = TE_id

    plt_rscript = '''

    TE_family_plotter <- function (te_dataframe, filename_id){
        rm(list=ls())

        library(tidyverse)
        library(patchwork)
        library(plotly)
        te_raw_data = te_dataframe
        output_id = filename_id
        
        te_raw_data$fpkm = as.numeric(te_raw_data$fpkm)
        te_raw_data$tpm = as.numeric(te_raw_data$tpm)
        te_raw_data$sample_tissue = as.factor(te_raw_data$sample_tissue)
        
        mean_tpm = mean(te_raw_data$tpm)
        max_tpm = max(te_raw_data$tpm)
        sample_number = nrow(te_raw_data)
          
        plt_TE_sample = ggplot(data=te_raw_data) +
          geom_segment(aes(x=sample_id, xend=sample_id, y=0, yend=tpm), color="#3366CC", alpha=0.4, size=0.8) +
          geom_point(aes(x=sample_id, y=tpm), color="#66CCCC", size=0.1) +
          geom_hline(yintercept = mean_tpm, color="red", linetype = "dashed", size=0.5) +
          geom_hline(yintercept = max_tpm, color="red", linetype = "dashed", size=0.5) +
          geom_text(aes(x=sample_number*2/3, y=mean_tpm*1.5), 
                    label=paste0("Mean TPM: ", round(mean_tpm,2))) +
          geom_text(aes(x=sample_number*2/3, y=max_tpm*0.9), 
                    label=paste0("Max TPM: ", round(max_tpm,2))) +
          theme_bw() +
          theme(
              panel.grid = element_blank(), 
        
              axis.line = element_line(colour = 'black'),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 10, hjust = 1),
              legend.position = 'none'
          ) +
          labs( x=paste0("samples (n=", sample_number, ")"), y="TPM")

        plt_TE_tissue = ggplot(data = te_raw_data, mapping = aes(x=sample_tissue, y=tpm)) + 
          geom_boxplot(alpha=0.5, fill="gray", outlier.fill="gray", outlier.size=0.3) +
          geom_jitter(color="#66CCCC", alpha=0.2, size=0.5, width=0.1) +
          theme_bw() +
          scale_fill_brewer() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position="none",
            panel.border = element_blank(),
          ) +
          labs(x = "Tissue types", y = "TPM")
          
        output_plt_TE_supfam <- plotly::ggplotly(plt_TE_sample)
        htmlwidgets::saveWidget(output_plt_TE_supfam, 
                                paste0("Mainapp/static/Temp_R_TE/", output_id,
                                       "_plt_TE_famsample.html"),
                                selfcontained = TRUE)
        
        output_plt_TE_group <- plotly::ggplotly(plt_TE_tissue)
        htmlwidgets::saveWidget(output_plt_TE_group, 
                                paste0("Mainapp/static/Temp_R_TE/", output_id,
                                       "_plt_TE_famtissue.html"),
                                selfcontained = TRUE)
    }

    '''

    robjects.r(plt_rscript)
    robjects.r['TE_family_plotter'](_te_df, TE_id)  # Execute function
    return TE_id
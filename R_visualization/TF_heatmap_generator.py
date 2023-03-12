import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import random

pandas2ri.activate()  # Used for variable conversion (pd.DataFrame)


def TF_tpm_heatmap(raw_tf_data):
    # Generate 10-digit random number for heatmap nomination
    output_filename = str(random.randint(1000000000, 9999999999)) + "heatmap.png"

    robjects.globalenv['tf_data'] = raw_tf_data
    robjects.globalenv['output_filename'] = output_filename

    plt_rscript_tpm = '''
    
    TF_tpm_heatmap <- function (raw_tf_data, output_filename){
        
        library(reshape2)  # data cleaning
        library(tidyverse)  # data visualization
        library(aplot)  # inserting barplot to the heatmap's side
        
        library(viridis) # Ornamental purposes
        library(ggthemes)
        library(RColorBrewer)
        
        tf_data$value = as.numeric(tf_data$tpm) # used for casting dataframe
        tf_data_cast = dcast(tf_data, gene_id~sample_id, fun.aggregate = sum)
        rownames(tf_data_cast) = tf_data_cast$gene_id
        tf_data_cast = tf_data_cast[,-1]
        
        
        tf_data$tpm = as.numeric(tf_data$tpm)
        
        heatmap = ggplot(data=tf_data, aes(x=sample_id, y=gene_id)) +
          geom_raster(alpha=0.5, aes(fill=tpm)) +
          
          theme_classic() +
          #scale_fill_gradient2_tableau() +
          scale_fill_gradient2(low="#003366", high="#0f1f1d", mid="#ffffff") +
          theme(axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank(),
                ) +
          geom_text(aes(label=signif(tpm,4), color=tpm), size=1.7)+
          scale_color_gradient2(low="#003366", high="black", mid="#dadad2") +
          guides(color=F) # hide legend for gradient text
        
        tf_bar = ggplot(tf_data,aes(x=1,y=gene_id))+
          geom_tile(aes(fill=tf_fam))+
          scale_x_continuous(expand = c(0,0))+
          theme(panel.background = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_text(angle = 0, size = 6),
                legend.position = "left",
                legend.title = element_blank()
                ) +
          scale_fill_brewer(palette = 2)
        
        sample_bar = ggplot(tf_data,aes(x=sample_id,y=1))+
          geom_tile(aes(fill=sample_tissue))+
          scale_x_continuous(expand = c(0,0))+
          theme(panel.background = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank(),
                axis.text.y = element_blank(),
                axis.text.x = element_text(angle = 15, size = 6),
                legend.position = "left",
                legend.title = element_blank()
                ) +
          scale_fill_brewer()
        
        pltmain = heatmap%>%
          insert_left(tf_bar,width = 0.05)%>%
          insert_bottom(sample_bar, height = 0.05)
        ggsave(paste0("Mainapp/static/Temp_R_heatmap/", output_filename), plot=pltmain, 
          width = ncol(tf_data_cast), height = nrow(tf_data_cast)/6, dpi = 300)
    }
    '''

    robjects.r(plt_rscript_tpm)
    robjects.r['TF_tpm_heatmap'](raw_tf_data, output_filename)  # Execute function

    return output_filename

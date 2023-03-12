import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

# File deletion
import shutil
from pathlib import Path

# Used for passing pandas dataframe variable to R
pandas2ri.activate()  # Used for variable conversion (pd.DataFrame)


def overview_barplot(exp_data, miRNA_id, query_format, query_species):
    # Detect and delete existing plots
    tmp_plot_path = "Mainapp/static/Temp_R_miRNA"
    for elm in Path(tmp_plot_path).glob(miRNA_id + '*'):
        elm.unlink() if elm.is_file() else shutil.rmtree(elm)  # delete folder and file

    robjects.globalenv['exp_data'] = exp_data
    plt_rscript_count = '''
        overview_barplot_T <- function (exp_data, miRNA_id){
          suppressMessages(library(ggplot2))
          suppressMessages(library(ggthemes))
          suppressMessages(library(plotly))
          suppressMessages(library(ggrepel))
          
          # Generate overview barplot using expression data from a single gene
          # sample_id Counts Normalized tissue_type
        
          exp_data$tissue_type <- as.factor(exp_data$tissue_type);
          exp_data$Counts <- as.numeric(as.character(exp_data$Counts));
          exp_data$Normalized <- as.numeric(as.character(exp_data$Normalized));
          exp_data <- exp_data %>% arrange(tissue_type, sample_id)
          # exp_data = exp_data %>% arrange(tissue_type, Counts)
          exp_data$id <- seq_len(nrow(exp_data))
        
        
          overview_barplot1 <- ggplot(data=exp_data) +
              geom_segment(aes(x=id, xend=id, y=0, yend=Counts),
                            color="#3366CC", alpha=0.4, size=1) +
              geom_point(aes(x=id, y=Counts, shape=sample_id), color="#66CCCC", size=1) +
              
              labs( x="Samples", y="Counts",caption = miRNA_id) +
              theme_bw() +
              theme(
                panel.border = element_blank(),
                legend.position="none",
              )
          output <- plotly::ggplotly(overview_barplot1)
          HTML_name <- paste0("Mainapp/static/Temp_R_miRNA", miRNA_id, "_plt1.html")
          htmlwidgets::saveWidget(output, HTML_name, selfcontained = TRUE)
          return(TRUE)
        }
        '''

    plt_rscript_normalized = '''
            overview_barplot_T <- function (exp_data, miRNA_id){
              suppressMessages(library(ggplot2))
              suppressMessages(library(ggthemes))
              suppressMessages(library(plotly))
              suppressMessages(library(ggrepel))

              # Generate overview barplot using expression data from a single gene
              # sample_id Counts Normalized tissue_type

              exp_data$tissue_type <- as.factor(exp_data$tissue_type);
              exp_data$Counts <- as.numeric(as.character(exp_data$Counts));
              exp_data$Normalized <- as.numeric(as.character(exp_data$Normalized));
              exp_data <- exp_data %>% arrange(tissue_type, sample_id)
              # exp_data = exp_data %>% arrange(tissue_type, Counts)
              exp_data$id <- seq_len(nrow(exp_data))


              overview_barplot1 <- ggplot(data=exp_data) +
                  geom_segment(aes(x=id, xend=id, y=0, yend=Normalized),
                                color="#3366CC", alpha=0.4, size=1) +
                  geom_point(aes(x=id, y=Normalized, shape=sample_id), color="#66CCCC", size=1) +

                  labs( x="Samples", y="Normalized",caption = miRNA_id) +
                  theme_bw() +
                  theme(
                    panel.border = element_blank(),
                    legend.position="none",
                  )
              output <- plotly::ggplotly(overview_barplot1)
              HTML_name <- paste0("Mainapp/static/Temp_R_miRNA", miRNA_id, "_plt1.html")
              htmlwidgets::saveWidget(output, HTML_name, selfcontained = TRUE)
              return(TRUE)
            }
            '''


    try:
        if query_format == 'Counts':
            if query_species in ["zea", "sorghum"]:
                robjects.r(plt_rscript_count)
                robjects.r['overview_barplot_T'](exp_data, miRNA_id)
            else:
                # Frontend query species error
                raise SystemError
        elif query_format == 'Normalized':
            if query_species in ["zea", "sorghum"]:
                robjects.r(plt_rscript_normalized)
                robjects.r['overview_barplot_T'](exp_data, miRNA_id)
            else:
                # Frontend query species error
                raise SystemError
        else:
            # Frontend query type error
            raise SystemError

    except SystemError:
        return None


def tissue_specific_barplot(exp_data, miRNA_id, query_format, query_species):

    # Existing plots already deleted in the above function, considering these two functions are used simultaneously.
    # If used individually, additional cleaning mechanisms should be added.

    robjects.globalenv['exp_data'] = exp_data
    plt2_rscript_count = '''
    
    tissue_specific_barplot_T <- function (exp_data, miRNA_id){
      
      suppressMessages(library(ggplot2))
      suppressMessages(library(ggthemes))
      suppressMessages(library(plotly))
      suppressMessages(library(ggrepel))
      
      # Generate tissue-specific barplot using expression data from a single gene
      # sample_id Counts tissue_type
    
      exp_data$tissue_type <- as.factor(exp_data$tissue_type);
      exp_data$Counts <- as.numeric(as.character(exp_data$Counts));
      exp_data <- exp_data %>% arrange(tissue_type, sample_id)
      # exp_data = exp_data %>% arrange(tissue_type, Counts)
      exp_data$id <- seq_len(nrow(exp_data))
    
    
      tissue_specific_barplot_t <- ggplot(data = exp_data, mapping = aes(x=tissue_type, y=Counts)) + 
            geom_boxplot(
                alpha=0.5,
                fill="gray",
                outlier.fill="gray",
                outlier.size=0.3
              ) +
              geom_jitter(
                color="#66CCCC",
                alpha=0.5,
                size=0.5,
                ) +
              theme_bw() +
              labs(
                x = "Tissue types",
                y = "Counts",
                caption = miRNA_id,
              ) +
              theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position="none",
                panel.border = element_blank(),
              )

    
      output <- plotly::ggplotly(tissue_specific_barplot_t)
      HTML_name <- paste0("Mainapp/static/Temp_R_miRNA/", miRNA_id, "_plt2.html")
      htmlwidgets::saveWidget(output, HTML_name, selfcontained = TRUE)
      return(TRUE)
    }
    '''
    plt2_rscript_normalized = '''

        tissue_specific_barplot_T <- function (exp_data, miRNA_id){

          suppressMessages(library(ggplot2))
          suppressMessages(library(ggthemes))
          suppressMessages(library(plotly))
          suppressMessages(library(ggrepel))

          # Generate tissue-specific barplot using expression data from a single gene
          # sample_id Counts Normalized tissue_type

          exp_data$tissue_type <- as.factor(exp_data$tissue_type);
          exp_data$Normalized <- as.numeric(as.character(exp_data$Normalized));
          exp_data <- exp_data %>% arrange(tissue_type, sample_id)
          # exp_data = exp_data %>% arrange(tissue_type, Normalized)
          exp_data$id <- seq_len(nrow(exp_data))


          tissue_specific_barplot_t <- ggplot(data = exp_data, mapping = aes(x=tissue_type, y=Normalized)) + 
                geom_boxplot(
                    alpha=0.5,
                    fill="gray",
                    outlier.fill="gray",
                    outlier.size=0.3
                  ) +
                  geom_jitter(
                    color="#66CCCC",
                    alpha=0.5,
                    size=0.5,
                    ) +
                  theme_bw() +
                  labs(
                    x = "Tissue types",
                    y = "Normalized",
                    caption = miRNA_id,
                  ) +
                  theme(
                    axis.text.x = element_text(angle = 45, hjust = 1),
                    legend.position="none",
                    panel.border = element_blank(),
                  )


          output <- plotly::ggplotly(tissue_specific_barplot_t)
          HTML_name <- paste0("Mainapp/static/Temp_R_miRNA/", miRNA_id, "_plt2.html")
          htmlwidgets::saveWidget(output, HTML_name, selfcontained = TRUE)
          return(TRUE)
        }
        '''

    try:
        if query_format == 'Counts':
            if query_species in ["zea", "sorghum"]:
                robjects.r(plt2_rscript_count)
                robjects.r['tissue_specific_barplot_T'](exp_data, miRNA_id)
            else:
                # Frontend query species error
                raise SystemError
        elif query_format == 'Normalized':
            if query_species in ["zea", "sorghum"]:
                robjects.r(plt2_rscript_normalized)
                robjects.r['tissue_specific_barplot_T'](exp_data, miRNA_id)
            else:
                # Frontend query species error
                raise SystemError

        else:
            # Frontend query type error
            return SystemError

    except SystemError:
        return None

#' Get data and summarize for main biotypes
#'
#' This function allows to parse given file, from XICRA biotype to produce a dataframe containing 
#' summary statistics for main RNA biotyoes (miRNA, tRNA, piRNA, miscRNA, lincRNA, etc)
#' @param data_biotypes Results file from XICRA biotype and/or featurecount + modifications
#' @keywords XICRA
#' @export
get_data <- function(data_biotypes) {
  library(reshape)
  library(stringr)

  ## read data
  # ---------------------------------------------------------------------------------
  ## miRNAseq
  biotypes <- read.csv(data_biotypes, sep=",", stringsAsFactors=T, header=T, row.names = 1)
  
  ## remove row
  biotypes.remove <- c("total")
  biotypes <- biotypes[!(row.names(biotypes) %in% biotypes.remove), ]
  
  ## types of RNA biotypes
  #rownames(biotypes)
  
  ## get new names
  biotypes$biotype<- ifelse(rownames(biotypes)=="miRNA", "miRNA",
                            ifelse(rownames(biotypes)=="piRNA", "piRNA", 
                                   ifelse(rownames(biotypes)=="tRNA", "tRNA", 
                                          ifelse(rownames(biotypes)=="lincRNA", "lincRNA", 
                                                 ifelse(rownames(biotypes)=="rRNA", "rRNA", 
                                                        ifelse(rownames(biotypes)=="processed_transcript", "Processed Transcript", 
                                                               ifelse(rownames(biotypes)=="protein_coding", "Protein Coding", 
                                                                      ifelse(rownames(biotypes)=="Unassigned_MultiMapping", "Align NoUniq", 
                                                                             ifelse(rownames(biotypes)=="multimapping", "Align NoUniq", 
                                                                                    ifelse(rownames(biotypes)=="Unassigned_Ambiguity", "No Feature", 
                                                                                           ifelse(rownames(biotypes)=="Unassigned_NoFeatures", "No Feature", 
                                                                                                  ifelse(rownames(biotypes)=="unmapped", "Not Align", 
                                                                                                         ifelse(rownames(biotypes)=="antisense", "Antisense", 
                                                                                                                ifelse(rownames(biotypes)=="misc_RNA", "misc_RNA", "Other"))))))))))))))
  ## replace NA
  biotypes[is.na(biotypes)] <- 0
  
  # ---------------------------------------------------------------------------------
  ## summarize
  # ---------------------------------------------------------------------------------
  sum.biotypes<- aggregate(biotypes[,-ncol(biotypes)], by=list(biotypes$biotype), "sum")
  sum.biotypes$Group.1 <- factor(sum.biotypes$Group.1, levels=unique(biotypes$biotype))
  
  ## generate percentage
  perc.biotypes<- cbind.data.frame(sum.biotypes$Group.1, 100 * prop.table( as.matrix(sum.biotypes[, 2:ncol(sum.biotypes)]), 2 ) )
  names(perc.biotypes)[1]<- "Biotypes"
  
  return(perc.biotypes)
}

#' Plot RNA biotype data into stack bar plots
#'
#' This function allows to plot given data from XICRA biotype and parse by get_data()
#' @param data_perc Data converted to 100% results.
#' @param palette_color_name Name of the color palette.
#' @param out_folder Output folder to store results.
#' @keywords XICRA
#' @export
plot_biotypes <- function(data_perc, palette_color_name="Set3", out_folder=NULL) {
 
  library(ggplot2)
  library(RColorBrewer)
  
  # ---------------------------------------------------------------------------------
  # Plot percentage Biotypes ggplot
  # ---------------------------------------------------------------------------------
  melt.biotypes <- melt(data_perc)
  melt.biotypes<- melt.biotypes[order(melt.biotypes$Biotypes, melt.biotypes$value,  decreasing=T),]
  names(melt.biotypes)[2] <- 'Samples'
  names(melt.biotypes)[3] <- 'Value'
  
  p <- ggplot(melt.biotypes, aes(x=Samples, y=Value, fill=Biotypes)) +
    geom_bar(stat="identity", colour="black") + 
    scale_fill_brewer(palette=palette_color_name) + theme_classic() + 
    theme(axis.text.x = element_text(face = "bold", size = 10, angle = 90))
  
  
  ## get message
  if (is.null(out_folder)){
    print("No printing plot. Provide a folder path")
  } else {
    
    ## save in pdf
    HCGB.IGTP.DAnalysis::save_pdf(folder_path = out_folder, name_file = "biotypes-plot", plot_given = p)
    
  }
  
  return(p)
}

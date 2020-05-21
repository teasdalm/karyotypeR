#! /usr/bin/env Rscript
# Rscript to Karyotype BAM files


## index the bam file if needed
index_bam <- function(input_file){
  index_file <- paste(input_file, ".bai", sep = "")
  if(file.exists(index_file)){
    message(index_file, " exits... skipping indexing")
  } else{
    message("Indexing ---> ", input_file)
    bai_command <- paste("samtools index", input_file)
    system(bai_command)
  }
}

# read idxstats
read_idxstats <- function(input_file){
  # get input files
  bam_file <- input_file
  bai_file <- paste0(input_file, ".bai")
  message("Reading ---> ", bam_file)

  # Do we need to index?
  index_bam(bam_file)

  # run idxstats command
  idx_command <- paste("samtools idxstats", bam_file)
  dat <-system(idx_command, intern=TRUE)

  # fmake output dataframe
  read_info <- data.frame(matrix
                          (unlist
                          (strsplit(dat, "\t")),
                          nrow=length(dat),
                          byrow=T),
                          stringsAsFactors=FALSE)

  # format nicely
  names(read_info) <- c("chrom", "length", "mapped_reads", "un_mapped_reads")
  read_info$length <- as.numeric(read_info$length)
  read_info$mapped_reads <- as.numeric(read_info$mapped_reads)
  read_info <- read_info[-nrow(read_info),]

  # return data
  return(read_info)
}

# plot data with ggplot2
plot_data <- function(idx_df, input_file = args[1]){
  base <- tools::file_path_sans_ext(input_file)
  output_file <- paste0(base, ".pdf")
  message("Plotting data to ---> ", output_file)
  library(ggplot2, quietly = TRUE)
  pdf(output_file)
  p <- ggplot(idx_df, aes(x=length, y=mapped_reads)) +
    geom_text(aes(label=chrom)) +
    xlab("Chromosome length (bp)") +
    ylab("Mapped reads") +
    geom_smooth(method=lm)
  print(p)
  garbage <- dev.off()
}

## calculate score as Skoglung et al (2015) wolf paper
calc_stat <- function(idx_df, chrom_1="chr1", chrom_x="chrX"){
  if(chrom_1 %in% idx_df$chrom & chrom_x %in% idx_df$chrom){
    chrom_1_stat <- idx_df[idx_df$chrom == chrom_1, "mapped_reads"] / idx_df[idx_df$chrom == chrom_1, "length"]
    chrom_x_stat <- idx_df[idx_df$chrom == chrom_x, "mapped_reads"] / idx_df[idx_df$chrom == chrom_x, "length"]
    skoglund_score <- chrom_x_stat/chrom_1_stat
    skoglund_score <- round(skoglund_score, 3)

    # have a guess about the sex
    if(skoglund_score >= 0.8){
      message("Skoglund et al (2015) score ---> ", skoglund_score, " ---> guess most likely female")
    } else if(skoglund_score <= 0.6 & skoglund_score >= 0.4){
      message("Skoglund et al (2015) score ---> ", skoglund_score, " ---> guess most likely male")
    } else{
      message("Skoglund et al (2015) score ---> ", skoglund_score, " ---> guess tricky look at graph")
    }
  } else {
    message("Can't find chromosome names skipping Skoglund score")
  }
}

## Entry point here

# Get command line args
args <- commandArgs(trailingOnly = TRUE)

# argument check
if(length(args) < 1 | length(args) > 3){
	message("##karyotpye.R##\nUsage Karyotype.R <file.bam> <Chromosome1 [chr1]> <ChromosomeX [chrX]>")
} else if (length(args) < 3 & length(args) > 1){
  message("please supply both none-default chromosome names for skoglund_score")
} else if (length(args) == 3){
  tab <- read_idxstats(args[1])
  plot_data(tab)
  calc_stat(tab, args[2], args[3])
} else if (length(args) == 1){
  tab <- read_idxstats(args[1])
  plot_data(tab)
  calc_stat(tab)
} else {
  message("##karyotpye.R##\nUsage Karyotype.R <file.bam> <Chromosome1 [chr1]> <ChromosomeX [chrX]>")
}



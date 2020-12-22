# ==========================================================================
# Libraries
# ==========================================================================
library('data.table')
library('ggplot2')
library('dplyr')
library('scales')
library('RColorBrewer')
library('argparser')
library("makeitprettier")


# ==========================================================================
# Arguments
# ==========================================================================
# Create a parser
p <- arg_parser("Plots the distribution of small RNA lengths and the first nucleotide composition.")

# Add command line arguments
p <- add_argument(p,
   "--tables",
   help="Results tables in the format output by summarizeNucleotideByReadLenght.py",
   type="character",
   nargs=Inf)

p <- add_argument(
   p,
   "--factors",
   help="path to file containing the normalization factors. It should be a tab-delimited file with two columns, one with the sample name and another with the normalization factor.",
   type="character",
   nargs=Inf)


# Parse the command line arguments
argv <- parse_args(p)
cov_files <- argv$tables
fac_path <- argv$factors


dir.create(file.path("figure"), showWarnings = FALSE)

# cov_files <- list.files(pattern='*.nuc_bias.txt')

read_and_summarize <- function(cov_file){
	message(print(cov_file))
	cov_tmp <- fread(cov_file)
	setnames(cov_tmp, c('Length', 'Nucleotide', 'Count'))
	cov_tmp[, Exp:=gsub('.nuc_bias.txt', '', cov_file), ]
   cov_tmp[, Strain:=gsub('^(\\w+)_\\w+_\\w+_\\w+', '\\1', Exp), ]
   cov_tmp[, Generation:=gsub('^\\w+_(\\w+)_\\w+_\\w+', '\\1', Exp), ]
   cov_tmp[, Treatment:=gsub('^\\w+_\\w+_(\\w+)_\\w+', '\\1', Exp), ]
   cov_tmp[, Replicate:=gsub('^\\w+_\\w+_\\w+_(\\w+)', '\\1', Exp), ]
	return(cov_tmp)
}

covs <- lapply(cov_files, read_and_summarize)
counts <- rbindlist(covs, fill=TRUE)

# fac_path <- '/fsimb/groups/imb-kettinggr/adomingues/projects/imb_ketting_2017_15_Ricardo_pid3/results/tracks/normalization_factors.txt'

message("Reading normalization factors")
fac <- fread(fac_path)
fac[, V1:=gsub('\\.\\w+$', '', V1),]
counts <- merge(counts, fac, by.x='Exp', by.y='V1')
counts[, Normalized_counts:=Count/V2 * 10^6, ]
counts[, Exp:=gsub("_RppH", "", Exp, ignore.case = TRUE), ]

ggplot(counts, aes(x=Length, y=Count, fill=Nucleotide)) +
   geom_bar(stat='identity') +
   scale_fill_prettier() +
   scale_y_continuous(labels=comma) +
   facet_grid(Strain ~ Replicate) + 
   theme_poster() +
   labs(
      x="Read length",
      y="Number of reads"
      )
ggsave("figure/nucleotide_bias_read_length.pdf")
ggsave("figure/nucleotide_bias_read_length.png")

ggplot(counts, aes(x=Length, y=Normalized_counts, fill=Nucleotide)) +
   geom_bar(stat='identity') +
   scale_fill_prettier() +
   scale_y_continuous(labels=comma) +
   facet_grid(Strain ~ Replicate) + 
   theme_poster() +
   labs(
      x="Read length",
      y="r.p.m.\n(mapped reads)"
      )
ggsave("figure/nucleotide_bias_read_length.normalized.pdf")
ggsave("figure/nucleotide_bias_read_length.normalized.png")

if ("Mapping" %in% colnames(counts)){
   ## only sense to use the space better
   ## 
   sense_counts <- counts[Mapping=="sense"]

   ggplot(sense_counts, aes(x=Length, y=Normalized_counts, fill=Nucleotide)) +
      geom_bar(stat='identity') +
      scale_fill_prettier() +
      scale_y_continuous(labels=comma) +
      facet_wrap(~ Exp, ncol=3) + 
      theme_poster() +
      labs(
         x="Read length",
         y="r.p.m.\n(mapped reads)"
         )
   ggsave("figure/nucleotide_bias_read_length.sense.normalized.pdf")
   ggsave("figure/nucleotide_bias_read_length.sense.normalized.png")


   # ==========================================================================
   # back to back plot
   # ==========================================================================

   commapos <- function (x, ...) 
   {
       format(abs(x), ..., big.mark = ",", scientific = FALSE, trim = TRUE)
   }

   # http://docs.ggplot2.org/current/geom_text.html
   df <- data.frame(
     x = c(max(counts$Length), max(counts$Length)),
     y = c(Inf, -Inf),
     text = c("Sense", "Antisense")
   )

   ggplot(counts, aes(x=Length)) +
      geom_bar(
         data=subset(counts, Mapping == "sense"),
         aes(y=Normalized_counts, fill=Nucleotide),
         stat = "identity") +
      geom_bar(
         data=subset(counts, Mapping == "antisense"),
         aes(y=-Normalized_counts, fill=Nucleotide),
         stat = "identity") +
      geom_text(data=df, aes(x=x, y=y, label = text), vjust = "inward", hjust = "inward") +
      facet_grid(Exp ~ .) + 
      scale_fill_prettier() +
      scale_y_continuous(labels=commapos) +
      geom_hline(yintercept=0, color="gray", linetype="dashed") + 
      theme_poster() +
      labs(
         x="Read length",
         y="r.p.m.\n(mapped reads)"
         ) +
      theme(
         axis.line.x = element_line(colour=NA),
         panel.border = element_rect(colour = 'gray'))
   ggsave("figure/nucleotide_bias_read_length.normalized.backtoback.pdf")
   ggsave("figure/nucleotide_bias_read_length.normalized.backtoback.png")



   # ==========================================================================
   # similar plot but with facets
   # ==========================================================================

   counts[,toPlot:=ifelse(Mapping == "antisense", -Normalized_counts, Normalized_counts)]
   counts[,Mapping:=factor(Mapping,
            levels = c("sense", "antisense")) ]

   ggplot(counts, aes(x=Length, y=toPlot, fill=Nucleotide)) +
      geom_bar(stat='identity') +
      scale_fill_prettier() +
      scale_y_continuous(labels=comma) +
      facet_grid(Mapping ~ Exp, scales = "free_y") + 
      theme_poster() +
      labs(
         x="Read length",
         y="r.p.m.\n(mapped reads)"
         )
   ggsave("figure/nucleotide_bias_read_length.normalized.backtoback_facet.pdf")
   ggsave("figure/nucleotide_bias_read_length.normalized.backtoback_facet.png")
}


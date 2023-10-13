library(ggplot2)
library(stringr)
library(argparse)
library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(tidyverse, warn.conflicts = FALSE)

#current_path = rstudioapi::getActiveDocumentContext()$path 
#setwd(dirname(current_path))
current_path = getwd()

# Parse arguments
parser <- ArgumentParser(description='Annotate Structural Variants using StructuralVariantAnnotation')
parser$add_argument("--vcf", type="character", required=TRUE,
                    help="VCF file from GRIDSS with structural variants")
parser$add_argument("--output_prefix", default=dirname(current_path), required=FALSE,
                    type="character",
                    help="Output path with prefix")
args <- parser$parse_args()

input_vcf <- file.path(args$vcf)
out_vcf <- file.path(paste(args$output_prefix, ".vcf", sep = ""))
out_bed <- file.path(paste(args$output_prefix, "_simple.bed", sep = ""))

#' Simple SV type classifier
simpleEventType <- function(gr) {
  pgr = partner(gr)
  return(ifelse(seqnames(gr) != seqnames(pgr), "CTX", # inter-chromosomosal
                ifelse(strand(gr) == strand(pgr), "INV",
                       ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
                              ifelse(xor(start(gr) < start(pgr), strand(gr) == "-"), "DEL",
                                     "DUP")))))
}

vcf <- VariantAnnotation::readVcf(input_vcf)
gr <- breakpointRanges(vcf)
gr <- gr[gr$FILTER == "PASS" & partner(gr)$FILTER == "PASS"]
svtype <- simpleEventType(gr)
info(vcf)$SIMPLE_TYPE <- NA_character_
info(vcf)$SVLEN <- NA_character_
info(vcf[gr$sourceId])$SIMPLE_TYPE <- svtype
info(vcf[gr$sourceId])$SVLEN <- gr$svLen
writeVcf(vcf, out_vcf)

simplegr <- gr[simpleEventType(gr) %in% c("INS", "INV", "DEL", "DUP")]
simplebed <- data.frame(
  chrom=seqnames(simplegr),
  # call the centre of the homology/inexact interval
  start=as.integer((start(simplegr) + end(simplegr)) / 2),
  end=as.integer((start(partner(simplegr)) + end(partner(simplegr))) / 2),
  name=simpleEventType(simplegr),
  score=simplegr$QUAL,
  length=simplegr$svLen,
  insSeq=simplegr$insSeq,
  strand="."
)
# Just the lower of the two breakends so we don't output everything twice
simplebed <- simplebed[simplebed$start < simplebed$end,]
write.table(simplebed, out_bed, quote=FALSE, sep='\t', 
            row.names=FALSE, col.names=FALSE)



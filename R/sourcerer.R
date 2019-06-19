#' Visualize chromosome bands
#'
#' @param path.to.cbmrscripts Absolute path to cbmrscripts
#' @return A series of functions
#' @export
source.bands <- function(path.to.cbmrscripts) {
    library(plotly)
    source(file.path(path.to.cbmrscripts, "scripts/bands/bands.R"), local=FALSE)
}
#' Extensive survival analysis for a series of signatures
#'
#' @param path.to.cbmrscripts Absolute path to cbmrscripts
#' @param path.to.tcgadump Absolute path to tcgadump
#' @param path.to.timer Absolute path to timer data
#' @return A series of functions
#' @export
source.survscan <- function(path.to.cbmrscripts,
                            path.to.tcgadump="/restricted/projectnb/montilab-p/personal/anthony/tcgadump",
                            path.to.timer="/restricted/projectnb/montilab-p/CBMrepositoryData/TCGA/tumorInfiltration/timer/TableS2.13059_2016_1028_MOESM3_ESM.txt") {
    library(Biobase)
    library(GSVA)
    library(ggplot2)
    library(gplots)
    library(dplyr)
    library(pheatmap)
    library(survival)
    library(survminer)
    source(file.path(path.to.cbmrscripts, "scripts/survscan/01_get_esets.R"), local=FALSE)
    source(file.path(path.to.cbmrscripts, "scripts/survscan/02_do_gsva.R"), local=FALSE)
    source(file.path(path.to.cbmrscripts, "scripts/survscan/03_immune_infiltration.R"), local=FALSE)
    source(file.path(path.to.cbmrscripts, "scripts/survscan/04_combine_signatures.R"), local=FALSE)
    source(file.path(path.to.cbmrscripts, "scripts/survscan/05_survival_analysis.R"), local=FALSE)
}
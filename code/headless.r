#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).txt", call. = FALSE)
}

source(args[1])
source("main.r")

warnings()

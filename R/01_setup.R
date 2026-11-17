# 01_setup.R
# Project-level configuration and path setup for the Ghana cocoa-market-access project.


project_root <- normalizePath(".")
raw_dir <- file.path(project_root, "data", "raw")
processed_dir <- file.path(project_root, "data", "processed")
outputs_dir <- file.path(project_root, "outputs")
figures_dir <- file.path(outputs_dir, "figures")
tables_dir <- file.path(outputs_dir, "tables")

if (!dir.exists(raw_dir)) stop("Raw data directory not found: ", raw_dir)
if (!dir.exists(processed_dir)) dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(readr)

options(scipen = 999)

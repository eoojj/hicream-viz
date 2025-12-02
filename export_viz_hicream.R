#! /bin/R

################################################################################
# Authors:
#   Sylvain Foissac, GenPhySE, Universite de Toulouse, INRAE, 31326 Castanet-Tolosan, France
#   Toby Dylan Hocking, LASSO lab, Département d'informatique, Universite de Sherbrooke, Sherbrooke, QC J1K-2R1, Canada
#   Elise Jorge, GenPhySE, Universite de Toulouse, INRAE, 31326 Castanet-Tolosan, France
#   Pierre Neuvial, Institut de Mathematiques de Toulouse, UMR 5219, CNRS UPS, 31062 Toulouse cedex 9, France
#   Nathalie Vialaneix, Universite de Toulouse, INRAE, MIAT, 31326 Castanet-Tolosan, France
#
# Copyright (C) 2022
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
################################################################################

# The purpose of this script is to produce dynamic visualization of hicream
# results for the three studied datasets and all studied resolutions and 
# chromosomes

## Load libraries
library(animint2)
library(data.table)
library(dplyr)
library(here)
library(IRanges)

## Load source functions
source("scripts/viz/source_viz_hicream.R")

## Options
options(scipen = 999)

# Create plots for all datasets, resolutions and chromosomes.

## Simu dataset
for(method in c("diff", "counts")){
  for(resolution in c(200000, 500000, 1000000)){
    for(chromosome in c(1, 7, 21)){
      res_path <- paste0("res/simu/", resolution, "/newexp_hicream_", method, "_chr", chromosome, "_", resolution, ".tsv")
      out_path <- paste0("res/viz/viz_hicream_simu_", method, "_chr", chromosome, "_", resolution)
      viz <- plot_dynamic_hicream(res_path)
      animint2dir(viz,  out.dir = out_path)
    }
  }
}

## Bonev dataset 
for(method in c("diff")){
  for(resolution in c(50000)){
    for(chromosome in c(1:19)){
      res_path <- paste0("res/bonev/", resolution, "/hicream_", method, "_chr", chromosome, "_", resolution, ".tsv")
      out_path <- paste0("res/viz/viz_hicream_bonev_", method, "_chr", chromosome, "_", resolution) 
      viz <- plot_dynamic_hicream(res_path)
      animint2dir(viz,  out.dir = out_path)
    }
  }
}

# CTCF dataset
for(method in c("diff")){
  for(resolution in c(100000, 200000)){
    # Create index of peaks from global index + ctcf peaks.
    index_path <- here(paste0("data/ctcf/", resolution, "/index.", resolution, ".bed4"))
    peaks_path <- here("data/ctcf/ctcf.peaks.mm9.bed")
    ctcf_path <- index_ctcf(index_path, peaks_path)
    for(chromosome in c(1:19)){
      res_path <- paste0("res/ctcf/", resolution, "/hicream_", method, "_chr", chromosome, "_", resolution, ".tsv")
      out_path <- paste0("res/viz/viz_hicream_ctcf_", method, "_chr", chromosome, "_", resolution)
      viz <- plot_dynamic_hicream(res_path, ctcf_path)
      animint2dir(viz,  out.dir = out_path)
    }
  }
}

sessionInfo()
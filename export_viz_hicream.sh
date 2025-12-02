#! /bin/bash

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

#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --mem=100G
#SBATCH -J export_viz_hicream
#SBATCH -o export_viz_hicream.log
#SBATCH -e Error_export_viz_hicream.log

#SBATCH --mail-type=BEGIN,END,FAIL 
#SBATCH --mail-user=elise.jorge@inrae.fr


module purge
module load statistics/R/4.4.0
module load tools/Pandoc/3.1.2
module load compilers/gcc/12.2.0

Rscript scripts/viz/export_viz_hicream.R

# Note: Run from the directory where the script is located with
# sbatch scripts/viz/export_viz_hicream.sh

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

# The purpose of this script is to create two functions: one to create index of
# ctcf peaks from Hi-C data index and raw ctcf peaks and one to generate dynamic
# animint2 visualization of hicream results

## Function to create index of ctcf peaks
index_ctcf <- function(index_path, peaks_path){
  index <- fread(here(index_path))
  names(index) <- c("chr", "start", "end", "bin")
  peaks <- fread(here(peaks_path))
  names(peaks) <- c("chr", "start", "end", "cond", "nb_peaks", ".")
  all_chr <- as.character(c(1:19))
  peaks <- peaks %>% filter(chr %in% all_chr)
  peaks$chr <- as.integer(peaks$chr)
  res <- lapply(c(1:19), function(chromosome){
    index_chr <- index %>% filter(chr==chromosome)
    peaks_chr <- peaks %>% filter(chr==chromosome)
    index_Iranges <- IRanges(start = index_chr$start, end = index_chr$end)
    peaks_Iranges <- IRanges(start = peaks_chr$start, end = peaks_chr$end)
    overlaps <- findOverlaps(index_Iranges, peaks_Iranges)
    res_chr <- data.frame(bin = index_chr[overlaps@from]$bin, nb_peaks =  peaks_chr[overlaps@to]$nb_peaks)
    res_chr <- res_chr %>% group_by(bin) %>% summarise(nb_peaks = sum(nb_peaks))
    res_chr$chr <- rep(chromosome, nrow(res_chr))
    return(res_chr[, c("chr", "bin", "nb_peaks")])
  })
  out_path <- paste0(strsplit(index_path, ".bed4")[[1]], ".ctcf_peaks.txt")
  write.table(do.call(rbind, res), out_path, row.names = FALSE)
  return(out_path)
}


## Function to create animint2 object from hicream results
plot_dynamic_hicream <- function(res_path, ctcf_path = NULL, pixel_cutoff = NULL){
  parts <- strsplit(basename(res_path), "_")[[1]] 
  chromosome <- gsub("chr", "", parts[grep("chr", parts)])
  resolution <- gsub(".tsv", "", parts[grep(".tsv", parts)])  
  pixel_dt <- fread(res_path)[, let(
    Cluster=factor(clust),
    neg.log10.p = -log10(p.value)
  )]
  if(!is.null(pixel_cutoff)){
    pixel_dt <- pixel_dt %>% filter(region1<=pixel_cutoff & region2<=pixel_cutoff)
  }
  
  if(!is.null(ctcf_path)){
    ctcf_index <- fread(ctcf_path, header = TRUE)[chr==chromosome]
    ##Normalize peaks to [0,10] for nice output
    ctcf_index$nb_peaks <-(ctcf_index$nb_peaks - min(ctcf_index$nb_peaks)) / (max(ctcf_index$nb_peaks) - min(ctcf_index$nb_peaks))
    ctcf_index$nb_peaks <- 10*ctcf_index$nb_peaks
    }
  
  r1r2_xy_mat <- rbind(
    c(0.5, -0.5),
    c(0.5, 0.5)) 
  set_xy <- function(DT, prefix){
    geti <- function(i)DT[[paste0(prefix,i)]]
    DT[, paste0(
      prefix, "_", c("x","y")
    ) := as.data.table(
      cbind(geti(1), geti(2)) %*% r1r2_xy_mat
    )][]
  }
  get_corners <- function(r1, r2, half.width){
    corners <- function(d1, d2)data.table(corner1=r1+d1, corner2=r2+d2)
    set_xy(rbind(
      corners(-half.width, -half.width),
      corners(-half.width, half.width),
      corners(half.width, half.width),
      corners(half.width, -half.width),
      corners(-half.width, -half.width)
    ), "corner")
  }
  expand <- 45
  expand.prop <- expand/100
  set_xy(pixel_dt, "region")
  pixel_corner_dt <- pixel_dt[, data.table(.SD, get_corners(region1, region2, expand.prop))]
  
  get_boundaries <- function(DT){
    dcast_input <- rbind(
      DT[, .(
        region1=seq(min(region1), max(region1)),
        region2=min(region2)-1L,
        value=0
      )],
      DT[, .(
        region1=min(region1)-1L,
        region2=seq(min(region2), max(region2)),
        value=0
      )],
      DT[, .(region1, region2, value=1)])
    wide <- dcast(dcast_input, region1 ~ region2, fill=0)
    m <- as.matrix(wide[,-1])
    clust_id_mat <- cbind(rbind(m, 0), 0)
    exp_one <- function(x)c(x, max(x)+1)
    path_list <- contourLines(
      exp_one(wide$region),
      exp_one(as.integer(colnames(m))),
      clust_id_mat, levels=0.5)
    xy <- c('x','y')
    circ_diff_vec <- function(z)diff(c(z,z[1]))
    circ_diff_dt <- function(DT, XY)DT[
      , paste0("d",XY) := lapply(.SD, circ_diff_vec), .SDcols=XY]
    out <- data.table(id=seq_along(path_list))[
      , path_list[[id]][xy]
      , by=id]
    for(XY in list(xy, paste0("d",xy))){
      circ_diff_dt(out, XY)
    }
    set_xy(out[
      , both_zero := ddx==0 & ddy==0
    ][!c(both_zero[.N], both_zero[-.N]), .(id, region1=x, region2=y)], "region")
  }
  
  myround <- function(x, bin_size=1, offset=0)round((x+offset)/bin_size)*bin_size
  off_list <- list(x=20, y=-24)
  for(xy in names(off_list)){
    round_fun <- function(rxy)myround(rxy, 50, off_list[[xy]])
    region_xy <- pixel_corner_dt[[paste0("region_",xy)]]
    corner_xy <- pixel_corner_dt[[paste0("corner_",xy)]]
    pixel_dt[, paste0("round_region_",xy) := round_fun(get(paste0("region_",xy)))]
    round_region_xy <- round_fun(region_xy)
    set(pixel_corner_dt, j=paste0("round_region_",xy), value=round_region_xy)
    set(pixel_corner_dt, j=paste0("rel_region_",xy), value=region_xy-round_region_xy)
    set(pixel_corner_dt, j=paste0("rel_corner_",xy), value=round(corner_xy-round_region_xy,2))
  }
  add_round_regions <- function(DT)DT[
    , round_regions := paste0(round_region_x,",",round_region_y)
  ][]
  add_round_regions(pixel_dt)
  add_round_regions(pixel_corner_dt)[, let(
    rel_regions = paste0(rel_region_x,",",rel_region_y)
  )][, .(corners=.N), by=rel_regions]
  
  local_cluster_presence <- unique(pixel_dt[, .(
    round_region_x, round_region_y, round_regions, Cluster)])
  local_cluster_dt <- pixel_dt[
    , get_boundaries(.SD)
    , by=.(Cluster,round_region_x,round_region_y,round_regions)]
  region_tiles <- pixel_corner_dt[, .(
    pixels=.N,
    mean.neg.log10.p=mean(neg.log10.p),
    mean.logFC=mean(logFC),
    clusters=length(unique(Cluster))
  ), by=.(round_regions,round_region_x,round_region_y)]
  for(xy in names(off_list)){
    local_cluster_dt[, paste0("rel_region_",xy) := get(paste0("region_",xy))-get(paste0("round_region_",xy))][]
    region_tiles[, (xy) := get(paste0("round_region_",xy))-off_list[[xy]]][]
  }
  
  corner_range <- dcast(
    pixel_corner_dt,
    . ~ .,
    list(min,max),
    value.var=c("rel_corner_x","rel_corner_y"))
  pixel_dt[, lnl.p := log10(neg.log10.p)][, let(
    round_logFC = myround(logFC,0.25),
    round_neg.log10.p = myround(neg.log10.p, 5),
    round_lnl.p = myround(lnl.p, 0.1)
  )][, let(
    volcano_bin = sprintf("logFC=%s -log10(p)=%s", round_logFC, round_neg.log10.p),
    lnl_bin = sprintf("logFC=%s l-l(p)=%s", round_logFC, round_lnl.p),
    relative_logFC = logFC-round_logFC,
    relative_neg.log10.p = neg.log10.p-round_neg.log10.p,
    relative_lnl.p = lnl.p-round_lnl.p
  )][]
  dcast(pixel_dt, . ~ ., list(min,max,Nvalues=function(x)length(unique(x))), value.var=patterns("round",cols=names(pixel_dt)))
  
  volcano_heat_lnl <- pixel_dt[, .(
    pixels=.N
  ), by=.(round_logFC,round_lnl.p,lnl_bin)]
  volcano_heat_lnl[order(pixels)]
  ggplot()+
    geom_tile(aes(
      round_logFC, round_lnl.p, fill=log10(pixels)),
      data=volcano_heat_lnl)+
    geom_text(aes(
      round_logFC, round_lnl.p, label=round(log10(pixels))),
      color="red",
      size=5,
      data=volcano_heat_lnl)+
    scale_fill_gradient(low="white",high="black")+
    theme_bw()
  
  volcano_heat <- pixel_dt[, .(
    pixels=.N
  ), by=.(round_logFC,round_neg.log10.p,volcano_bin)]
  volcano_heat[order(pixels)]
  ggplot()+
    geom_tile(aes(
      round_logFC, round_neg.log10.p, fill=log10(pixels)),
      data=volcano_heat)+
    geom_text(aes(
      round_logFC, round_neg.log10.p, label=round(log10(pixels))),
      color="red",
      size=5,
      data=volcano_heat)+
    scale_fill_gradient(low="white",high="black")+
    theme_bw()
  
  ## TODO break up huge log10(p)=0 counts into negative space bins.
  
  count_by_Cluster <- dcast(
    pixel_dt,
    Cluster ~ .,
    list(min, max),
    value.var=c("logFC", "neg.log10.p")
  )[
    pixel_dt[, .(
      pixels=.N,
      tiles=length(unique(round_regions))
    ), keyby=Cluster]
  ]
  count_by_Cluster_tile <- pixel_dt[, .(
    displayed_pixels=.N
  ), keyby=.(Cluster, round_regions)][
    count_by_Cluster
  ][, let(
    label=sprintf(
      "Cluster %s: %d/%d pixels shown, logFC %.1f to %.1f, log10(p) %.1f to %.1f",
      Cluster, displayed_pixels, pixels,
      logFC_min, logFC_max,
      neg.log10.p_min, neg.log10.p_max),
    rel_x = corner_range[, (rel_corner_x_max+rel_corner_x_min)/2],
    rel_y = corner_range$rel_corner_y_max+1
  )]
  first_list <- as.list(pixel_dt[which.max(neg.log10.p), .(volcano_bin, Cluster)])
  first_list$round_regions <- local_cluster_dt[Cluster==first_list$Cluster, round_regions][1]
  
  table(count_by_Cluster$tiles)
  hist(log10(count_by_Cluster$pixels))
  
  cluster_heat_dt <- count_by_Cluster[, let(
    log2.tiles=log2(tiles),
    log10.pixels=log10(pixels)
  )][, let(
    round.log2.tiles=myround(log2.tiles, 0.25),
    round.log10.pixels=myround(log10.pixels, 0.25)
  )][, let(
    size_bin=sprintf(
      "log2(tiles)=%.2f log10(pixels)=%.2f",
      round.log2.tiles, round.log10.pixels),
    rel.log2.tiles=log2.tiles-round.log2.tiles,
    rel.log10.pixels=log10.pixels-round.log10.pixels
  )][, .(
    clusters=.N
  ), by=.(size_bin, round.log2.tiles, round.log10.pixels)]
  ggplot()+
    geom_tile(aes(
      round.log10.pixels, round.log2.tiles,
      fill=log10(clusters)),
      data=cluster_heat_dt)+
    scale_fill_gradient(low="white",high="black")+
    theme_bw()
  
  count_by_Cluster[
    , rel_row := 1:.N
    , by=.(size_bin, rel.log10.pixels)]
  ggplot()+
    geom_point(aes(
      rel.log10.pixels, rel_row),
      data=count_by_Cluster)+
    facet_wrap("size_bin")
  
  ###Create CTCF data 
  if(!is.null(ctcf_path)){
    pixel_corner_dt_diag <- pixel_corner_dt %>% filter(region1==region2)
    comp_ctcf_data <- lapply(1:nrow(region_tiles), function(i){
      ##For the region, vertical projection of pixels 
      sel_region <- region_tiles[i,] %>% select(round_regions)
      zone <- paste0(strsplit(sel_region$round_regions, ",")[[1]][1],",0") ##vertical projection
      ##Select diagonal pixels of the projection
      sel_pixels <- pixel_corner_dt_diag %>% filter(round_regions == zone)
      ##Data frame with bin numbers and x start/end coordinates
      sel_pixels <- sel_pixels %>%
        group_by(region1) %>%
        summarise(
          min_x = min(rel_corner_x),
          max_x = max(rel_corner_x)
        )
      ##Add nb_peaks for each bin
      ctcf_data <- sel_pixels 
      ctcf_data$round_regions <- rep(sel_region$round_regions, nrow(ctcf_data))
      names(ctcf_data)[1] <- c("bin")
      ctcf_data$nb_peaks <- ctcf_index[match(ctcf_data$bin, ctcf_index$bin),]$nb_peaks
      return(ungroup(ctcf_data))
    })
    
    comp_ctcf_data <- bind_rows(comp_ctcf_data)
    
  }
  
  cluster.color <- "green"
  viz.common <- animint(
    title= paste0("Visualization of hicream results for chr", chromosome, " at ", as.integer(resolution)/1000, "kb resolution"),
    first=first_list,
    volcanoHeat=ggplot()+
      ggtitle("Volcano plot summary")+
      geom_tile(aes(
        round_logFC, round_neg.log10.p, fill=log10(pixels)),
        color="grey",
        size=1,
        data=volcano_heat)+
      geom_point(aes(
        logFC, neg.log10.p),
        data=pixel_dt,
        color="green",
        fill=cluster.color,
        chunk_vars="Cluster",
        showSelected="Cluster")+
      geom_tile(aes(
        round_logFC, round_neg.log10.p),
        clickSelects="volcano_bin",
        fill="transparent",
        data=volcano_heat)+
      scale_fill_gradient(low="white",high="black")+
      theme_bw()+
      theme_animint(width=300, height=300),
    genomeZoom=ggplot()+
      ggtitle("Genomic interaction zoom")+
      geom_polygon(aes(
        rel_corner_x, rel_corner_y,
        fill=TPRate,
        group=rel_regions),
        chunk_vars = "round_regions", #add option when files are too large
        clickSelects="Cluster",
        alpha_off=1,
        alpha=1,
        data=pixel_corner_dt,
        showSelected="round_regions", 
        color="grey") +
      {
        if (!is.null(ctcf_path)) {
          geom_rect(aes( xmin = min_x,
                         xmax = max_x, 
                         ymin = -nb_peaks,
                         ymax = 0), 
                    fill = "grey60", 
                    color = "black", 
                    alpha = 0.4, 
                    data = comp_ctcf_data, 
                    showSelected = "round_regions")
        } else {
          NULL
        }
      } +
      geom_path(aes(
        rel_region_x, rel_region_y, group=paste(Cluster, id)),
        color=cluster.color,
        alpha=1,
        alpha_off=0.2,
        clickSelects="Cluster",
        showSelected="round_regions",
        data=local_cluster_dt)+
      geom_text(aes(
        rel_x, rel_y, label=label, group=1),
        size=20,
        data=count_by_Cluster_tile,
        chunk_vars="round_regions",
        showSelected=c("Cluster","round_regions"))+
      scale_fill_gradient2(
        low = "white",
        high = "deeppink3", 
        name = "TDP"
      ) + 
      theme_bw()+
      theme_animint(height=800, width=800, rowspan=3),
    clusterHeat=ggplot()+
      ggtitle("Cluster size summary")+
      geom_tile(aes(
        round.log10.pixels, round.log2.tiles,
        fill=log10(clusters)),
        color="grey",
        data=cluster_heat_dt)+
      geom_point(aes(
        log10.pixels, log2.tiles),
        color="black",
        fill=cluster.color,
        showSelected="Cluster",
        data=count_by_Cluster)+
      geom_tile(aes(
        round.log10.pixels, round.log2.tiles),
        fill="transparent",
        color="black",
        clickSelects="size_bin",
        data=cluster_heat_dt)+
      scale_fill_gradient(low="white",high="black")+
      theme_bw()+
      theme_animint(width=300, height=300, last_in_row=TRUE),
    volcanoZoom=ggplot()+
      ggtitle("Volcano plot zoom")+
      geom_point(aes(
        relative_logFC, relative_neg.log10.p),
        data=pixel_dt,
        size=4,
        showSelected="volcano_bin",
        chunk_vars="volcano_bin",
        fill="white",
        color=cluster.color,
        color_off="black",
        clickSelects="Cluster")+
      theme_bw()+
      theme_animint(width=300, height=300),
    clusterZoom=ggplot()+
      ggtitle("Cluster size zoom")+
      geom_point(aes(
        rel.log10.pixels, rel_row),
        size=4,
        fill="white",
        color=cluster.color,
        color_off="black",
        showSelected="size_bin",
        clickSelects="Cluster",
        data=count_by_Cluster)+
      scale_fill_gradient(low="white",high="black")+
      theme_bw()+
      theme_animint(width=300, height=300, last_in_row=TRUE),
    genomeSummary=ggplot()+
      ggtitle("Genomic interaction summary")+
      geom_tile(aes(
        x, y, fill=mean.logFC),
        data=region_tiles,
        color="transparent")+
      geom_point(aes(
        round_region_x-off_list$x,
        round_region_y-off_list$y,
        group=round_regions),
        data=local_cluster_presence,
        color="black",
        fill=cluster.color,
        showSelected="Cluster")+
      geom_tile(aes(
        x, y),
        data=region_tiles,
        clickSelects="round_regions",
        fill="transparent",
        color="black")+
      scale_x_continuous(
        "Genomic bin")+
      scale_y_continuous(
        "Interaction distance")+
      scale_fill_gradient2()+
      theme_bw()+
      theme_animint(width=300, height=300), 
    source = "https://github.com/eoojj/hicream-viz/export_viz_hicream.R"
    
  )
  return(viz.common)
}

sessionInfo()
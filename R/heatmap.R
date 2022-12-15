#' Build heatmap matrix
#'
#' @param cl_mat edge-list for cluster to visualise
#' @param res resolution of the heatmap
#' @param val_col column with colormap values
#'
#' @return Sparse Matrix with cluster values for heatmap visualisation
#' @export
#'
full_f_mat<-function(cl_mat,res,val_col){

  range_5kb<-range(unique(c(cl_mat$ego,cl_mat$alter)))
  bin_5kb<-seq(range_5kb[1],range_5kb[2],by=res)
  #add the bins not present in original Hi-C dataset
  #miss_bin<-bin_5kb[which(!(bin_5kb %in% unique(c(mat_df$X1,mat_df$X2))))]

  id_conv<-seq_along(bin_5kb)
  names(id_conv)<-bin_5kb

  cl_mat$ego_id<-id_conv[as.character(cl_mat$ego)]
  cl_mat$alter_id<-id_conv[as.character(cl_mat$alter)]

  #chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=sqrt(-log10(cl_mat$pois.pval)),symmetric = T)
  chr_mat<-Matrix::sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=cl_mat[[val_col]],symmetric = T,dimnames = list(names(id_conv),names(id_conv)))
  return(chr_mat)
}

#' Utility function to produce color maps
#'
#' @param res_set Set of resolution over which to map colors
#' @param palette Palette to use for color map
#'
#' @return list with color map hexdecimal values and corresponding breaks
#' @export
#'
custom_color_map_fn<-function(res_set,palette){
  p_breaks<-c(unlist(lapply(seq_along(res_set),function(x){
    seq((x-1)*100,(x-1)*100 +100,length.out = 101)[-101]
  })),(length(res_set)-1)*100 +100)

  p_color<-rev(RColorBrewer::brewer.pal(n=length(res_set),name = palette))
  p_col<-unlist(lapply(seq_along(res_set),function(x){
    grDevices::colorRampPalette(c("black",p_color[x]))(100)
  }))
  return(list(p.col=p_col,p.breaks=p_breaks))
}

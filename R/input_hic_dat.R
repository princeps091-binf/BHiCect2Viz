
get_mres_hic_dat<-function(dat_folder,chromo){
  chr_dat_l<-lapply(hub_res_set,function(x)read_delim(file = paste0(dat_folder,x,'/',chromo,'.txt'),delim = '\t',col_names = F))
  names(chr_dat_l)<-hub_res_set
  chr_dat_l<-lapply(chr_dat_l,function(x){
    tmp<-x%>%
      filter(!(is.nan(X3)))%>%
      filter(X1!=X2)
    tmp_bin<-unique(c(x$X1,x$X2))
    tmp_self<-tibble(X1=tmp_bin,X2=tmp_bin) %>%
      left_join(.,x)
    min_self<-min(tmp_self$X3,na.rm=T)
    tmp_self<-tmp_self %>%
      mutate(X3=ifelse(is.na(X3),min_self,X3))
    return(tmp %>%
             bind_rows(.,tmp_self))
  })
  chr_dat_l<-lapply(chr_dat_l,function(x){
    preprocessParams <- BoxCoxTrans(x$X3,na.rm = T)
    x <- data.frame(x,weight=predict(preprocessParams, x$X3))
    x$weight<-x$weight+(1-min(x$weight,na.rm = T))
    return(x)
  })

}
# Produce color-scale separating each resolution into separate color-channel

make_mres_color_map<-function(hub_res_set,res_set,chr_dat_l){
  out_chr_dat_l<-lapply(seq_along(hub_res_set),function(x){
    tmp_dat<-chr_dat_l[[hub_res_set[x]]]
    res_idx<-which(res_set == hub_res_set[x])
    toMin<-(res_idx-1)*100 +1
    toMax<-(res_idx-1)*100 +99
    tmp_dat$color<-toMin+(tmp_dat$weight-min(tmp_dat$weight))/(max(tmp_dat$weight)-min(tmp_dat$weight))*(toMax-toMin)
    return(tmp_dat)
  })
  names(out_chr_dat_l)<-hub_res_set
  return(out_chr_dat_l)
}

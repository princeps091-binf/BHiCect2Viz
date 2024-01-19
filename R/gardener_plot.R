# Plotgardener routine
#' Heatmap plotting with plotgardener
#'
#' @param cl_col_dat Tibble collecting all heatmap entries to plot with corresponding coordinates
#' @param col_specs List collecting the necessary components to build color-map
#' @param chromo String indicating the chromosome containing the considered cluster
#' @param genome String indicating the genome symbol containing the considered cluster
#' @param figure_height Height in inches of the heatmap to plot
#' @param figure_width Width in inches of the heatmap to plot
#' @param res_num Named vector with the resolutions and their labels present in the HiC data
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @return Heatmap plot produced by the plotgardener backend
#' @export
#'
heatmap_plotgardener_fn<-function(cl_col_dat,col_specs,chromo,genome, figure_height, figure_width,res_num){


  plotgardener::pageCreate(width = figure_width, height =figure_height, default.units = "inches",showGuides = F)
  page_viewport<-grid::current.viewport()

  test<-plotgardener::plotHicSquare(
    resolution = 5000,
    data = cl_col_dat,
    chrom = chromo,
    chromstart = min(c(cl_col_dat$ego,cl_col_dat$alter)),
    chromend = max(c(cl_col_dat$ego,cl_col_dat$alter)),
    assembly = genome,
    x = figure_width/20, y = figure_height/20, width = figure_width/2, height = figure_height/2, default.units = "inches",
    draw=F)

  heat_grob_id<-names(test[["grobs"]][["children"]])[2]
  tmp_ego<-test[["grobs"]][["children"]][[heat_grob_id]]$x
  tmp_alter<-test[["grobs"]][["children"]][[heat_grob_id]]$y
  tmp_col<-test[["grobs"]][["children"]][[heat_grob_id]]$gp$fill

  tmp_hic_col<-tibble::tibble(ego=as.numeric(stringr::str_split_fixed(tmp_ego,"native",n=2)[,1]),
                      alter=as.numeric(stringr::str_split_fixed(tmp_alter,"native",n=2)[,1]),
                      fill=tmp_col)

  tmp_col<-col_specs$p.col

  names(tmp_col)<-as.character(1:length(tmp_col))

  tmp_conv<-cl_col_dat %>%
    dplyr::bind_rows(tibble::tibble(ego=cl_col_dat$alter,
                            alter=cl_col_dat$ego,
                            color=cl_col_dat$color)) %>%
    dplyr::mutate(code=tmp_col[as.character(round(.data$color))]) %>%
    dplyr::distinct()

  tmp_hic_col<-tmp_hic_col %>%
    dplyr::left_join(tmp_conv)

  test[["grobs"]][["children"]][[heat_grob_id]]$gp$fill<-tmp_hic_col$code

  # Plotting
  plotgardener::pagePlotPlace(
    plot=test,
    x = figure_width/4, y = figure_width/4,
    width = figure_width/2, height = figure_height/2,default.units = "inches"
  )

  ideogramPlot <- plotgardener::plotIdeogram(
    chrom = chromo, assembly = genome,
    showBands=F,
    orientation = "h",
    x = figure_width/4, y = figure_width/4 - figure_width/10,
    width = figure_width/2, height = figure_height/25, just = "left"
  )

  region <- plotgardener::pgParams(chrom = chromo,
                     chromstart = min(c(cl_col_dat$ego,cl_col_dat$alter)),
                     chromend = max(c(cl_col_dat$ego,cl_col_dat$alter)))
  plotgardener::annoHighlight(
    plot = ideogramPlot, params = region,
    fill = "red",
    y = figure_width/4 - figure_width/10, height = figure_height/20, just = c("centre", "center"), default.units = "inches"
  )
  plotgardener::annoZoomLines(
    plot = ideogramPlot, params = region,
    y0 = figure_width/4 - figure_width/10 + figure_height/25,
    x1 = c(figure_width/4, figure_width/4 + figure_width/2),
    y1 = figure_height/4, default.units = "inches"
  )
  plotgardener::annoGenomeLabel(
    plot = test,
    axis = "y",
    x = figure_width/4 - figure_width/20, y = figure_width/4,
    width = figure_width/2, height = figure_width/20,
    default.units = "inches",scale = "Mb"
  )
  grid::seekViewport(page_viewport$name)
  col_scale<- matrix(rep(rev(col_specs$p.col),50),ncol=50)
  grid::grid.raster(col_scale,
                    x=grid::unit(figure_width/4 + figure_width/2 + figure_width/25  ,"inches"),
                    y=grid::unit(figure_width/2,"inches"),
                    height=grid::unit(0.5,"npc"),
                    width=grid::unit(0.05,"npc"))

  for(i in seq_along(res_num)){
    grid::grid.text(
      label = paste(names(res_num)[i]),
      x=grid::unit(figure_width/4 + figure_width/2 + figure_width/15  ,"inches"),
      y=grid::unit(figure_width/4 + i*figure_width/12,"inches"),
      just = c("left", "top")
    )

  }

}

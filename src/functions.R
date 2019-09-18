write_xlsx <- function(df, excel_file, sheet, overwritesheet=TRUE, ...){
  df <- as.data.frame(df)
  if (!file.exists(excel_file)) {
    xlsx::write.xlsx(df, excel_file, sheet, row.names=FALSE, ...)
    return(as_tibble(df))
  }
  wb <- xlsx::loadWorkbook(excel_file)  
  sheets <- readxl::excel_sheets(excel_file)
  if (sheet %in% sheets){
    if (overwritesheet==FALSE) stop(base::simpleError(glue::glue("Sheet {sheet} exists in file and 'overwritesheet' set to FALSE, please use xlsx::write.xlsx() with a new sheet name and 'append=TRUE' option")))
    xlsx::removeSheet(wb, sheet)
    xlsx::saveWorkbook(wb, excel_file)
    xlsx::write.xlsx(df, excel_file, sheet, row.names=FALSE, ...)
    return(as_tibble(df))
  }
  xlsx::write.xlsx(df, excel_file, sheet, row.names=FALSE, ...)
  return(as_tibble(df))
}


givemedonutsorgivemedeath <- function(dataset, group, subgroup, values, legend=TRUE, legend_offset=c(-0.1,-0.15), inner_labs=TRUE, inner_labs_col="white", cols=c('firebrick2','forestgreen','dodgerblue2'), outer_cols="shades") {
  #   ## house keeping
  #   if (missing(file)) file <- getwd()
  #   plot.new(); op <- par(no.readonly = TRUE); on.exit(par(op))
  #
  #   pdf(file, width = width, height = height, bg = 'snow')
  #
  
  #   dataset=family_summary
  #   values="log_count"
  #   group="phylum"
  #   subgroup="family"
  ## useful values and colors to work with
  ## each group will have a specific color
  ## each subgroup will have a specific shade of that color
  nr <- nrow(dataset)
  width <- max(sqrt(dataset[[values]])) / 0.8
  #tbl <- with(browsers, table(browser)[order(unique(browser))])
  tbl <- table(dataset[[group]])[unique(dataset[[group]])]
  cols <- unlist(Map(rep, cols, tbl))
  alphas <- list()
  for (i in names(tbl)) {
    alphas[[i]] <- seq.int(from=200, by=(-15), length.out =tbl[[i]])
  }
  
  ## loop creates pie slices
  plot.new()
  par(omi = c(0.5,0.5,0.75,0.5), mai = c(0.1,0.1,0.1,0.1), las = 1)
  for (i in 1:nr) {
    par(new = TRUE)
    
    ## create color/shades
    rgb <- col2rgb(cols[i])
    f0 <- rep(NA, nr)
    f0[i] <- rgb(rgb[1], rgb[2], rgb[3], unlist(alphas)[i], maxColorValue = 255)
    
    ## stick labels on the outermost section
    lab <- dataset[[subgroup]]
    if (dataset[[values]][i] == max(dataset[[values]])) {
      lab0 <- lab
    } else lab0 <- NA
    
    ## plot the outside pie and shades of subgroups
    pie(dataset[[values]], border = NA, radius = 5 / width, col = f0,
        labels = lab0, cex = 1.25)
    
    ## repeat above for the main groups
    par(new = TRUE)
    rgb <- col2rgb(cols[i])
    f0[i] <- rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
    
    pie(dataset[[values]], border = NA, radius = 4 / width, col = f0, labels = NA)
  }
  
  if(legend){
    # par(mar=c(5.2, 4.1, 4.1, 8.2), xpd=TRUE)
    #legend_offset=c(-0.2,-0.15)
    legend("topright",inset=legend_offset,legend=unique(dataset[[group]]),pch=rep(15,'.',length(dataset[[group]])), pt.cex=1.2, pt.lwd=2, border=unique(cols), col=unique(cols),bty='n', title = paste0(toupper(substr(group, 1, 1)), substr(group, 2, nchar(group))), title.adj = 0.3,cex=1.1)
    #text("topright", labels = paste0(toupper(substr(group, 1, 1)), substr(group, 2, nchar(group))), adj=0.3,cex=1.2)
  }
  # pch=rep(15,'.',length(pctr)),
  ## extra labels on graph
  if (inner_labs){
    sums <- aggregate(dataset$count_num, list(dataset[[group]]), sum) %>% arrange(desc(x))
    text(x = c(-.05, .25, 0.73), y = c(.15, -.40, -.09),
         labels = paste0("n = ",prettyNum(sums[,2], big.mark = ",", preserve.width="none")), col = inner_labs_col, cex = 1.25)
  }
}

plot_horiz_facet <- function(dataset, group.by, facet.by){
  
}
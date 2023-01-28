#' matriannotate: a function to add matrisome annotations to gene tables
#'
#' @param data character. The input data.
#' @param gene.column character. The name of the input column with gene names (accepted formats are Gene Symbols, NCBI gene ID (formerly Entrez Gene ID) and Ensembl Gene ID).
#' @param species character. One of the species for which matrisome annotations are provided (human, mouse, c.elegans, drosophila, zebrafish).
#'
#' @return The input data with added matrisome categories and families.
#'
#' @examples matriannotate(data, gene.column, species)

matriannotate <- function(data=NULL,
                          gene.column=NULL,
                          species=NULL){
  if(is.null(data)){
    cat(crayon::red("no data provided, execution stops \n"))
  }
  if(class(data)!="data.frame"){
    cat(crayon::red("data should be in data.frame format, execution stops \n"))
  }
  if(is.null(gene.column)){
    cat(crayon::red("a column indicating gene IDs must be provided, execution stops \n"))
  }
  if(is.null(species)){
    cat(crayon::red("no species provided, execution stops \n"))
  }

  df <- data
  n <- gene.column
  if(species=="human"){k <- matrisome.list$human}
  if(species=="mouse"){k <- matrisome.list$mouse}
  if(species=="c.elegans"){k <- matrisome.list$c.elegans
  k$family <- ifelse(k$family=="ECM-affiliated","ECM-affiliated Proteins",as.character(k$family))
  }
  if(species=="zebrafish"){k <- matrisome.list$zebrafish
  k$category <- ifelse(k$category=="not.available","Non-matrisome",as.character(k$category))
  k$family <- ifelse(k$family=="","Non-matrisome",
                     ifelse(k$family=="not.available","Non-matrisome",as.character(k$family)))
  }
  if(species=="drosophila"){k <- matrisome.list$drosophila
  k$category <- ifelse(k$category=="Homologs/Orthologs to Mammalian Matrisome-Associated Genes","Matrisome-associated",
                       ifelse(k$category=="Homologs/Orthologs to Mammalian Core Matrisome Genes","Core matrisome","Drosophila matrisome"))
  k$family <- ifelse(k$family=="ECM-affiliated","ECM-affiliated Proteins",as.character(k$family))
  }

  df2 <- distinct(merge(k,df,by.x="gene",by.y=as.character(n),all.y = T))
  names(df2)[1:3] <- c("Annotated Gene",
                       "Annotated Matrisome Division",
                       "Annotated Matrisome Category")
  df2$`Annotated Gene`[df2$`Annotated Gene`==""] <- "gene name missing in original data"
  df2$`Annotated Matrisome Division`[is.na(df2$`Annotated Matrisome Division`)] <- "Non-matrisome"
  df2$`Annotated Matrisome Category`[is.na(df2$`Annotated Matrisome Category`)] <- "Non-matrisome"
  df2[is.na(df2)] <- ""
  attributes(df2)$workflow <- "matrisomeannotatoR"
  return(df2)
}

#' matrianalyze: tabulations of matriannotated data
#'
#' @param data character. The matriannotated input data.
#'
#' @return a tabulated data.frame.
#'
#' @examples matriannotate(data)

matrianalyze <- function(data=NULL){
  if(is.null(data)){
    cat(crayon::red("no data provided, execution stops \n"))
  }
  if(class(data)!="data.frame"){
    cat(crayon::red("data should be in data.frame format, execution stops \n"))
  }
  if(length(attributes(data))<1){
    cat(crayon::red("data should be annotated first, execution stops \n"))
  }
  if(attributes(data)$workflow!="matrisomeannotatoR"){
    cat(crayon::red("data should be annotated first, execution stops \n"))
  }

  n <- "Annotated Gene"
  df <- data
  tr <- suppressWarnings(
    apply(df,2,function(x){
      v <- gsub("[^0-9.-]", "",x)
      as.numeric(v)
    })
  )
  colnames(tr) <- colnames(df)
  tr <- tr[,!(is.na(colSums(tr)))]
  nmtr <- colnames(tr)
  tr <- as.data.frame(cbind(df[,colnames(df)%in%as.character(n)],tr))
  tr[,c(2:ncol(tr))] <- apply(tr[,c(2:ncol(tr))],2,as.numeric)
  colnames(tr) <- c(as.character(n),nmtr)
  df2_2 <- df[,c(1:3)]

  bf <- distinct(merge(df2_2,tr,by.x="Annotated Gene",by.y=as.character(n)))
  a <- aggregate(bf[4:ncol(bf)],list(bf$`Annotated Matrisome Division`),sum)
  b <- aggregate(bf[4:ncol(bf)],list(bf$`Annotated Matrisome Category`),sum)
  z <- bind_rows(a,b)
  names(z)[1] <- "Matrisome Annotation"
  attributes(z)$workflow <- "matrisomeanalyzeR"

  return(z)
}

#' matribar: barplots for matriannotated data
#'
#' @param data character. The matriannotated input data.
#' @param print.plot logical. If TRUE (default) the plot is printed, otherwise the graph object (in its original format) is returned.
#'
#' @return a ggplot2 object (printed or not)
#'
#' @examples matribar(data)

matribar <- function(data=NULL,
                     print.plot=TRUE){
  if(is.null(data)){
    cat(crayon::red("no data provided, execution stops \n"))
  }
  if(class(data)!="data.frame"){
    cat(crayon::red("data should be in data.frame format, execution stops \n"))
  }
  if(attributes(data)$workflow!="matrisomeannotatoR"){
    cat(crayon::red("graphs can only be drawn for annotated files, execution stops \n"))
  }

  gdf <- data
  v1 <- as.data.frame(table(gdf$`Annotated Matrisome Division`))
  v1$source <- "Annotated Matrisome Division"
  v2 <- as.data.frame(table(gdf$`Annotated Matrisome Category`))
  v2$source <- "Annotated Matrisome Category"
  d1 <- bind_rows(v1,v2)

  d1$Var1 <- factor(d1$Var1,levels = c("Drosophila matrisome",
                                       "Nematode-specific core matrisome",
                                       "Nematode-specific matrisome-associated",
                                       "Putative Matrisome",
                                       "Core matrisome",
                                       "Matrisome-associated",
                                       "Apical Matrix",
                                       "Cuticular Collagens",
                                       "Cuticlins",
                                       "ECM Glycoproteins",
                                       "Collagens",
                                       "Proteoglycans",
                                       "ECM-affiliated Proteins",
                                       "ECM Regulators",
                                       "Secreted Factors",
                                       "Non-matrisome"
  ))
  nd1 <- data.frame(V1=c(
    "Drosophila matrisome",
    "Nematode-specific core matrisome",
    "Nematode-specific matrisome-associated",
    "Putative Matrisome",
    "Core matrisome",
    "Matrisome-associated",
    "Apical Matrix",
    "Cuticular Collagens",
    "Cuticlins",
    "ECM Glycoproteins",
    "Collagens",
    "Proteoglycans",
    "ECM-affiliated Proteins",
    "ECM Regulators",
    "Secreted Factors",
    "Non-matrisome"),
    V2=c("#B118DB",
         "#B118DB",
         "#741B47",
         "#B118DB",
         "#002253",
         "#DB3E18",
         "#B3C5FC",
         "#B3C5FC",
         "#BF9000",
         "#13349D",
         "#0584B7",
         "#59D8E6",
         "#F4651E",
         "#F9A287",
         "#FFE188",
         "#D9D9D9"
    )
  )
  nd1 <- nd1[nd1$V1%in%d1$Var1,2]

  p1 <- ggplot(d1,aes(Var1,Freq,fill=Var1)) +
    geom_bar(stat="Identity") +
    scale_fill_manual(values=as.character(nd1)) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    xlab("") + ylab("Counts") +
    facet_wrap(~d1$source,scales="free")

  if(isTRUE(print.plot)){
  suppressMessages(
    suppressWarnings(
      plot(p1)
    )
  )}else{
    return(p1)
  }
}

#' matriflow: alluvials for matriannotated data
#'
#' @param data character. The matriannotated input data.
#' @param print.plot logical. If TRUE (default) the plot is printed, otherwise the graph object (in its original format) is returned.
#'
#' @return a ggplot2 object (printed or not)
#'
#' @examples matriflow(data)

matriflow <- function(data=NULL,
                      print.plot=TRUE){
  if(is.null(data)){
    cat(crayon::red("no data provided, execution stops \n"))
  }
  if(class(data)!="data.frame"){
    cat(crayon::red("data should be in data.frame format, execution stops \n"))
  }
  if(attributes(data)$workflow!="matrisomeannotatoR"){
    cat(crayon::red("graphs can only be drawn for annotated files, execution stops \n"))
  }

  require(ggalluvial)
  gdf <- data
  d2 <- as.data.frame(table(gdf$`Annotated Matrisome Division`,gdf$`Annotated Matrisome Category`))
  d2 <- d2[d2$Freq>0,]
  d2[is.na(d2)] <- "Non-matrisome"
  d2$Var1 <- factor(d2$Var1,levels=c("Drosophila matrisome",
                                     "Nematode-specific core matrisome",
                                     "Nematode-specific matrisome-associated",
                                     "Putative Matrisome",
                                     "Core matrisome",
                                     "Matrisome-associated",
                                     "Non-matrisome"))
  d2$Var2 <- factor(d2$Var2,levels=c("Apical Matrix",
                                     "Cuticular Collagens",
                                     "Cuticlins",
                                     "ECM Glycoproteins",
                                     "Collagens",
                                     "Proteoglycans",
                                     "ECM-affiliated Proteins",
                                     "ECM Regulators",
                                     "Secreted Factors",
                                     "Non-matrisome"))
  vc <- data.frame(all=c("Drosophila matrisome",
          "Nematode-specific core matrisome",
          "Nematode-specific matrisome-associated",
          "Putative Matrisome",
          "Core matrisome",
          "Matrisome-associated",
          "Apical Matrix",
          "Cuticular Collagens",
          "Cuticlins",
          "ECM Glycoproteins",
          "Collagens",
          "Proteoglycans",
          "ECM-affiliated Proteins",
          "ECM Regulators",
          "Secreted Factors",
          "Non-matrisome"
          ),
          colz=c(
            "#B118DB",
            "#B118DB",
            "#741B47",
            "#B118DB",
            "#002253",
            "#DB3E18",
            "#B3C5FC",
            "#B3C5FC",
            "#BF9000",
            "#13349D",
            "#0584B7",
            "#59D8E6",
            "#F4651E",
            "#F9A287",
            "#FFE188",
            "#D9D9D9"
          ))
  vc <- vc[vc$all%in%unique(c(d2$Var1,d2$Var2)),]
  vc <- vc[vc$all%in%d2$Var2,2]

  p2 <- ggplot(d2,aes(y = Freq, axis1 = Var1, axis2 = Var2)) +
    geom_alluvium(aes(fill = Var2), width = 1/12) +
    geom_stratum(width = 1/12, fill = "white", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Annotated Matrisome Division", "Annotated Matrisome Category"), expand = c(.05, .05)) +
    scale_fill_manual(values = as.character(vc)) +
    theme_bw() + theme(legend.position = "none",
                       axis.text.y = element_blank(),
                       axis.ticks.y = element_blank()) + ylab("")
  if(isTRUE(print.plot)){
  suppressMessages(
    suppressWarnings(
      plot(p2)
    )
  )}else{
    return(p2)
  }
}

#' matripie: donutpie for matriannotated data
#'
#' @param data character. The matriannotated input data.
#'
#' @return a DonutPie object (always printed to screen)
#'
#' @examples matripie(data)

matripie <- function(data=NULL){
  if(is.null(data)){
    cat(crayon::red("no data provided, execution stops \n"))
  }
  if(class(data)!="data.frame"){
    cat(crayon::red("data should be in data.frame format, execution stops \n"))
  }
  if(attributes(data)$workflow!="matrisomeannotatoR"){
    cat(crayon::red("graphs can only be drawn for annotated files, execution stops \n"))
  }
  require(moonBook)

  gdf <- data
  d2 <- as.data.frame(table(gdf$`Annotated Matrisome Division`,gdf$`Annotated Matrisome Category`))
  d2 <- d2[d2$Freq>0,]
  d2$Var1 <- factor(d2$Var1,levels=c("Drosophila matrisome",
                                     "Nematode-specific core matrisome",
                                     "Nematode-specific matrisome-associated",
                                     "Putative Matrisome",
                                     "Core matrisome",
                                     "Matrisome-associated",
                                     "Non-matrisome"))
  d2$Var2 <- factor(d2$Var2,levels=c("Apical Matrix",
                                     "Cuticular Collagens",
                                     "Cuticlins",
                                     "ECM Glycoproteins",
                                     "Collagens",
                                     "Proteoglycans",
                                     "ECM-affiliated Proteins",
                                     "ECM Regulators",
                                     "Secreted Factors",
                                     "Non-matrisome"))

  PieDonutCustom <- function(data, mapping, start = getOption("PieDonut.start",
                                                                   0), addPieLabel = TRUE, addDonutLabel = TRUE, showRatioDonut = TRUE,
                                  showRatioPie = TRUE, ratioByGroup = TRUE, showRatioThreshold = getOption("PieDonut.showRatioThreshold",
                                                                                                           0.02), labelposition = getOption("PieDonut.labelposition",
                                                                                                                                            2), labelpositionThreshold = 0.1, r0 = getOption("PieDonut.r0",
                                                                                                                                                                                             0.3), r1 = getOption("PieDonut.r1", 1), r2 = getOption("PieDonut.r2",
                                                                                                                                                                                                                                                    1.2), explode = NULL, selected = NULL, explodePos = 0.1,
                                  color = "white", pieAlpha = 0.8, donutAlpha = 1, maxx = NULL,
                                  showPieName = TRUE, showDonutName = FALSE, title = NULL,
                                  pieLabelSize = 4, donutLabelSize = 3, titlesize = 5, explodePie = TRUE,
                                  explodeDonut = FALSE, use.label = TRUE, use.labels = TRUE,
                                  family = getOption("PieDonut.family", ""), palette_name="Dark2")
      {
        (cols = colnames(data))
        if (use.labels)
          data = moonBook::addLabelDf(data, mapping)
        count <- NULL
        if ("count" %in% names(mapping))
          count <- moonBook::getMapping(mapping, "count")
        count
        pies <- donuts <- NULL
        (pies = moonBook::getMapping(mapping, "pies"))
        if (is.null(pies))
          (pies = moonBook::getMapping(mapping, "pie"))
        if (is.null(pies))
          (pies = moonBook::getMapping(mapping, "x"))
        (donuts = moonBook::getMapping(mapping, "donuts"))
        if (is.null(donuts))
          (donuts = moonBook::getMapping(mapping, "donut"))
        if (is.null(donuts))
          (donuts = moonBook::getMapping(mapping, "y"))
        if (!is.null(count)) {
          df <- data %>% group_by(.data[[pies]]) %>% dplyr::summarize(Freq = sum(.data[[count]]))
          df
        }
        else {
          df = data.frame(table(data[[pies]]))
        }
        colnames(df)[1] = pies
        df$end = cumsum(df$Freq)
        df$start = dplyr::lag(df$end)
        df$start[1] = 0
        total = sum(df$Freq)
        df$start1 = df$start * 2 * pi/total
        df$end1 = df$end * 2 * pi/total
        df$start1 = df$start1 + start
        df$end1 = df$end1 + start
        df$focus = 0
        if (explodePie)
          df$focus[explode] = explodePos
        df$mid = (df$start1 + df$end1)/2
        df$x = ifelse(df$focus == 0, 0, df$focus * sin(df$mid))
        df$y = ifelse(df$focus == 0, 0, df$focus * cos(df$mid))
        df$label = df[[pies]]
        df$ratio = df$Freq/sum(df$Freq)
        if (showRatioPie) {
          df$label = ifelse(df$ratio >= showRatioThreshold, paste0(df$label,
                                                                   "\n(", scales::percent(df$ratio), ")"),
                            as.character(df$label))
        }
        df$labelx = (r0 + r1)/2 * sin(df$mid) + df$x
        df$labely = (r0 + r1)/2 * cos(df$mid) + df$y
        if (!is.factor(df[[pies]]))
          df[[pies]] <- factor(df[[pies]])
        df
        mainCol = RColorBrewer::brewer.pal(nrow(df), name=palette_name)
        df$radius = r1
        df$radius[df$focus != 0] = df$radius[df$focus != 0] + df$focus[df$focus !=
                                                                         0]
        df$hjust = ifelse((df$mid%%(2 * pi)) > pi, 1, 0)
        df$vjust = ifelse(((df$mid%%(2 * pi)) < (pi/2)) | (df$mid%%(2 *
                                                                      pi) > (pi * 3/2)), 0, 1)
        df$segx = df$radius * sin(df$mid)
        df$segy = df$radius * cos(df$mid)
        df$segxend = (df$radius + 0.05) * sin(df$mid)
        df$segyend = (df$radius + 0.05) * cos(df$mid)
        df
        if (!is.null(donuts)) {
          subColor = makeSubColor(mainCol, no = length(unique(data[[donuts]])))
          subColor
          data
          if (!is.null(count)) {
            df3 <- as.data.frame(data[c(donuts, pies, count)])
            colnames(df3) = c("donut", "pie", "Freq")
            df3
            df3 <- eval(parse(text = "complete(df3,donut,pie)"))
            df3$Freq[is.na(df3$Freq)] = 0
            if (!is.factor(df3[[1]]))
              df3[[1]] = factor(df3[[1]])
            if (!is.factor(df3[[2]]))
              df3[[2]] = factor(df3[[2]])
            df3 <- df3 %>% arrange(.data$pie, .data$donut)
            a <- df3 %>% spread(.data$pie, value = .data$Freq)
            a = as.data.frame(a)
            a
            rownames(a) = a[[1]]
            a = a[-1]
            a
            colnames(df3)[1:2] = c(donuts, pies)
          }
          else {
            df3 = data.frame(table(data[[donuts]], data[[pies]]),
                             stringsAsFactors = FALSE)
            colnames(df3)[1:2] = c(donuts, pies)
            a = table(data[[donuts]], data[[pies]])
            a
          }
          a
          df3
          df3$group = rep(colSums(a), each = nrow(a))
          df3$pie = rep(1:ncol(a), each = nrow(a))
          total = sum(df3$Freq)
          total
          df3$ratio1 = df3$Freq/total
          df3
          if (ratioByGroup) {
            df3$ratio = scales::percent(df3$Freq/df3$group)
          }
          else {
            df3$ratio <- scales::percent(df3$ratio1)
          }
          df3$end = cumsum(df3$Freq)
          df3
          df3$start = dplyr::lag(df3$end)
          df3$start[1] = 0
          df3$start1 = df3$start * 2 * pi/total
          df3$end1 = df3$end * 2 * pi/total
          df3$start1 = df3$start1 + start
          df3$end1 = df3$end1 + start
          df3$mid = (df3$start1 + df3$end1)/2
          df3$focus = 0
          if (!is.null(selected)) {
            df3$focus[selected] = explodePos
          }
          else if (!is.null(explode)) {
            selected = c()
            for (i in 1:length(explode)) {
              start = 1 + nrow(a) * (explode[i] - 1)
              selected = c(selected, start:(start + nrow(a) -
                                              1))
            }
            selected
            df3$focus[selected] = explodePos
          }
          df3
          df3$x = 0
          df3$y = 0
          df
          if (!is.null(explode)) {
            explode
            for (i in 1:length(explode)) {
              xpos = df$focus[explode[i]] * sin(df$mid[explode[i]])
              ypos = df$focus[explode[i]] * cos(df$mid[explode[i]])
              df3$x[df3$pie == explode[i]] = xpos
              df3$y[df3$pie == explode[i]] = ypos
            }
          }
          df3$no = 1:nrow(df3)
          df3$label = df3[[donuts]]
          if (showRatioDonut) {
            if (max(nchar(levels(df3$label))) <= 2)
              df3$label = paste0(df3$label, "(", df3$ratio,
                                 ")")
            else df3$label = paste0(df3$label, "\n(", df3$ratio,
                                    ")")
          }
          df3$label[df3$ratio1 == 0] = ""
          df3$label[df3$ratio1 < showRatioThreshold] = ""
          df3$hjust = ifelse((df3$mid%%(2 * pi)) > pi, 1, 0)
          df3$vjust = ifelse(((df3$mid%%(2 * pi)) < (pi/2)) | (df3$mid%%(2 *
                                                                           pi) > (pi * 3/2)), 0, 1)
          df3$no = factor(df3$no)
          df3
          labelposition
          if (labelposition > 0) {
            df3$radius = r2
            if (explodeDonut)
              df3$radius[df3$focus != 0] = df3$radius[df3$focus !=
                                                        0] + df3$focus[df3$focus != 0]
            df3$segx = df3$radius * sin(df3$mid) + df3$x
            df3$segy = df3$radius * cos(df3$mid) + df3$y
            df3$segxend = (df3$radius + 0.05) * sin(df3$mid) +
              df3$x
            df3$segyend = (df3$radius + 0.05) * cos(df3$mid) +
              df3$y
            if (labelposition == 2)
              df3$radius = (r1 + r2)/2
            df3$labelx = (df3$radius) * sin(df3$mid) + df3$x
            df3$labely = (df3$radius) * cos(df3$mid) + df3$y
          }
          else {
            df3$radius = (r1 + r2)/2
            if (explodeDonut)
              df3$radius[df3$focus != 0] = df3$radius[df3$focus !=
                                                        0] + df3$focus[df3$focus != 0]
            df3$labelx = df3$radius * sin(df3$mid) + df3$x
            df3$labely = df3$radius * cos(df3$mid) + df3$y
          }
          df3$segx[df3$ratio1 == 0] = 0
          df3$segxend[df3$ratio1 == 0] = 0
          df3$segy[df3$ratio1 == 0] = 0
          df3$segyend[df3$ratio1 == 0] = 0
          if (labelposition == 0) {
            df3$segx[df3$ratio1 < showRatioThreshold] = 0
            df3$segxend[df3$ratio1 < showRatioThreshold] = 0
            df3$segy[df3$ratio1 < showRatioThreshold] = 0
            df3$segyend[df3$ratio1 < showRatioThreshold] = 0
          }
          df3
          del = which(df3$Freq == 0)
          del
          if (length(del) > 0)
            subColor <- subColor[-del]
          subColor
        }
        p <- ggplot() + ggforce::theme_no_axes() + coord_fixed()
        if (is.null(maxx)) {
          r3 = r2 + 0.3
        }
        else {
          r3 = maxx
        }
        p1 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x", y0 = "y",
                                                   r0 = as.character(r0), r = as.character(r1), start = "start1",
                                                   end = "end1", fill = pies), alpha = pieAlpha, color = color,
                                        data = df) + transparent() + scale_fill_manual(values = mainCol) +
          xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
        if ((labelposition == 1) & (is.null(donuts))) {
          p1 <- p1 + geom_segment(aes_string(x = "segx",
                                             y = "segy", xend = "segxend", yend = "segyend"),
                                  data = df) + geom_text(aes_string(x = "segxend",
                                                                    y = "segyend", label = "label", hjust = "hjust",
                                                                    vjust = "vjust"), size = pieLabelSize, data = df,
                                                         family = family)
        }
        else if ((labelposition == 2) & (is.null(donuts))) {
          p1 <- p1 + geom_segment(aes_string(x = "segx",
                                             y = "segy", xend = "segxend", yend = "segyend"),
                                  data = df[df$ratio < labelpositionThreshold, ]) +
            geom_text(aes_string(x = "segxend", y = "segyend",
                                 label = "label", hjust = "hjust",
                                 vjust = "vjust"), size = pieLabelSize,
                      data = df[df$ratio < labelpositionThreshold,
                      ], family = family) + geom_text(aes_string(x = "labelx",
                                                                 y = "labely", label = "label"), size = pieLabelSize,
                                                      data = df[df$ratio >= labelpositionThreshold, ],
                                                      family = family)
        }
        else {
          p1 <- p1 + geom_text(aes_string(x = "labelx", y = "labely",
                                          label = "label"), size = pieLabelSize, data = df,
                               family = family)
        }
        if (showPieName)
          p1 <- p1 + annotate("text", x = 0, y = 0, label = pies,
                              size = titlesize, family = family)
        p1 <- p1 + theme(text = element_text(family = family))
        if (!is.null(donuts)) {
          if (explodeDonut) {
            p3 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x",
                                                       y0 = "y", r0 = as.character(r1), r = as.character(r2),
                                                       start = "start1", end = "end1", fill = "no",
                                                       explode = "focus"), alpha = donutAlpha,
                                            color = color, data = df3)
          }
          else {
            p3 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x",
                                                       y0 = "y", r0 = as.character(r1), r = as.character(r2),
                                                       start = "start1", end = "end1", fill = "no"),
                                            alpha = donutAlpha, color = color, data = df3)
          }
          p3 <- p3 + transparent() + scale_fill_manual(values = subColor) +
            xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
          p3
          if (labelposition == 1) {
            p3 <- p3 + geom_segment(aes_string(x = "segx",
                                               y = "segy", xend = "segxend", yend = "segyend"),
                                    data = df3) + geom_text(aes_string(x = "segxend",
                                                                       y = "segyend", label = "label", hjust = "hjust",
                                                                       vjust = "vjust"), size = donutLabelSize,
                                                            data = df3, family = family)
          }
          else if (labelposition == 0) {
            p3 <- p3 + geom_text(aes_string(x = "labelx",
                                            y = "labely", label = "label"), size = donutLabelSize,
                                 data = df3, family = family)
          }
          else {
            p3 <- p3 + geom_segment(aes_string(x = "segx",
                                               y = "segy", xend = "segxend", yend = "segyend"),
                                    data = df3[df3$ratio1 < labelpositionThreshold,
                                    ]) + geom_text(aes_string(x = "segxend",
                                                              y = "segyend", label = "label", hjust = "hjust",
                                                              vjust = "vjust"), size = donutLabelSize,
                                                   data = df3[df3$ratio1 < labelpositionThreshold,
                                                   ], family = family) + geom_text(aes_string(x = "labelx",
                                                                                              y = "labely", label = "label"), size = donutLabelSize,
                                                                                   data = df3[df3$ratio1 >= labelpositionThreshold,
                                                                                   ], family = family)
          }
          if (!is.null(title))
            p3 <- p3 + annotate("text", x = 0, y = r3,
                                label = title, size = titlesize, family = family)
          else if (showDonutName)
            p3 <- p3 + annotate("text", x = (-1) * r3,
                                y = r3, label = donuts, hjust = 0, size = titlesize,
                                family = family)
          p3 <- p3 + theme(text = element_text(family = family))
          grid::grid.newpage()
          print(p1, vp = grid::viewport(height = 1, width = 1))
          print(p3, vp = grid::viewport(height = 1, width = 1))
        }
        else {
          p1
        }
      }
      d2$Var1 <- as.character(d2$Var1)
      d2$Var2 <- as.character(d2$Var2)
      suppressMessages(
        suppressWarnings(
      PieDonutCustom(d2,aes(Var1,Var2,Freq),palette_name = "Dark2")
    )
  )
}

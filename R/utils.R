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

  suppressMessages(require(dplyr))

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
    stop()
  }
  if(class(data)!="data.frame"){
    cat(crayon::red("data should be in data.frame format, execution stops \n"))
    stop()
  }
  if(attributes(data)$workflow!="matrisomeannotatoR"){
    cat(crayon::red("graphs can only be drawn for annotated files, execution stops \n"))
    stop()
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
    stop()
  }
  if(class(data)!="data.frame"){
    cat(crayon::red("data should be in data.frame format, execution stops \n"))
    stop()
  }
  if(attributes(data)$workflow!="matrisomeannotatoR"){
    cat(crayon::red("graphs can only be drawn for annotated files, execution stops \n"))
    stop()
  }

  suppressMessages(require(ggalluvial))
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
  names(d2)[2] <- "Annotated Matrisome Categories"

  p2 <- ggplot(d2,aes(y = Freq, axis1 = Var1, axis2 = `Annotated Matrisome Categories`)) +
    geom_alluvium(aes(fill = `Annotated Matrisome Categories`), width = 1/12) +
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

#' matriring: donut plot for matriannotated data
#'
#' @param data character. The matriannotated input data.
#' @param print.plot logical. If TRUE (default) the plot is printed, otherwise the graph object (in its original format) is returned.
#'
#' @return a ggplot2 object (printed or not)
#'
#' @examples matriring(data)

matriring <- function(data=NULL,
                      print.plot=TRUE){

  if(is.null(data)){
    cat(crayon::red("no data provided, execution stops \n"))
    stop()
  }
  if(class(data)!="data.frame"){
    cat(crayon::red("data should be in data.frame format, execution stops \n"))
    stop()
  }
  if(attributes(data)$workflow!="matrisomeannotatoR"){
    cat(crayon::red("graphs can only be drawn for annotated files, execution stops \n"))
    stop()
  }

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
  fin <- distinct(merge(d2,vc,by.x="Var2",by.y="all"))
  fin$Var1 <- NULL
  fin$Var2 <- factor(fin$Var2,levels=c("Apical Matrix",
                                       "Cuticular Collagens",
                                       "Cuticlins",
                                       "ECM Glycoproteins",
                                       "Collagens",
                                       "Proteoglycans",
                                       "ECM-affiliated Proteins",
                                       "ECM Regulators",
                                       "Secreted Factors",
                                       "Non-matrisome"))
  fin$lab <- paste0(fin$Var2," (",fin$Freq,")")
  fin <- fin[order(fin$Var2),]
  fin$ymax <- cumsum(fin$Freq)
  fin$ymin <- c(0, head(fin$ymax, n=-1))

  zc <- fin
  zc <- zc[order(zc$Var2),]

  p <- suppressWarnings(suppressMessages(
    ggplot(fin, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var2)) +
      geom_rect() +
      scale_fill_manual(values = as.character(zc$colz),
                        name="Annotated Matrisome Categories",
                        labels=zc$lab) +
      coord_polar(theta="y") +
      xlim(c(-1, 4)) +
      theme_void()
  ))

  if(isTRUE(print.plot)){
    suppressMessages(
      suppressWarnings(
        plot(p)
      )
    )}else{
      return(p)
    }
}

#' matristar: polar bars for matriannotated data
#'
#' @param data character. The matriannotated input data.
#' @param print.plot logical. If TRUE (default) the plot is printed, otherwise the graph object (in its original format) is returned.
#'
#' @return a ggplot2 object (printed or not)
#'
#' @examples matristar(data)

matristar <- function(data=NULL,
                      print.plot=TRUE){

  if(is.null(data)){
    cat(crayon::red("no data provided, execution stops \n"))
    stop()
  }
  if(class(data)!="data.frame"){
    cat(crayon::red("data should be in data.frame format, execution stops \n"))
    stop()
  }
  if(attributes(data)$workflow!="matrisomeannotatoR"){
    cat(crayon::red("graphs can only be drawn for annotated files, execution stops \n"))
    stop()
  }

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
  fin <- distinct(merge(d2,vc,by.x="Var2",by.y="all"))
  fin$Var1 <- as.factor(as.character(fin$Var1))
  fin$labs <- paste0(fin$Var2," (",fin$Freq,")")

  if(fin$Freq[fin$Var1=="Non-matrisome"] > 4*max(fin$Freq[fin$Var1!="Non-matrisome"]) ){
    cat(crayon::red("warning: the Non-matrisome bar was scaled down to aid visualization"))
    fin$Freq[fin$Var1=="Non-matrisome"] <- 2*max(fin$Freq[fin$Var1!="Non-matrisome"])
  }else{
    fin<-fin
  }

  empty_bar <- 1
  to_add <- data.frame( matrix(NA, empty_bar*nlevels(fin$Var1), ncol(fin)) )
  colnames(to_add) <- colnames(fin)
  to_add$Var1 <- rep(levels(fin$Var1), each=empty_bar)
  data <- rbind(fin, to_add)
  data <- data %>% arrange(Var1)
  data$id <- seq(1, nrow(data))

  label_data <- data
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar
  label_data$hjust <- ifelse( angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)

  base_data <- data %>%
    group_by(Var1) %>%
    summarize(start=min(id), end=max(id) - empty_bar) %>%
    rowwise() %>%
    mutate(title=mean(c(start, end)))
  base_data <- base_data[base_data$Var1!="Non-matrisome",]

  p <- suppressMessages(suppressWarnings(
    ggplot(data, aes(x=as.factor(id), y=Freq, fill=colz)) +
      geom_bar(stat="identity",color="black") +
      ylim(-100,120) +
      theme_minimal() +
      scale_fill_manual(values = as.character(data$colz)) +
      theme(
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-1,4), "cm")
      ) +
      coord_polar() +
      geom_text(data=label_data, aes(x=id, y=Freq+10, label=labs, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
      geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )
  ))

  if(isTRUE(print.plot)){
    suppressMessages(
      suppressWarnings(
        plot(p)
      )
    )}else{
      return(p)
    }
}


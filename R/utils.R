#' matriannotate: a function to add matrisome annotations to gene tables
#'
#' @param data character. The input data.
#' @param gene.column character. The name of the input column with gene names (accepted formats are Gene Symbols, NCBI gene ID (formerly Entrez Gene ID) and Ensembl Gene ID).
#' @param species character. One of the species for which matrisome annotations are provided (human, mouse, c.elegans, drosophila, zebrafish, quail).
#'
#' @return The input data with added matrisome categories and families.
#'
#' @examples matriannotate(data, gene.column, species)

matriannotate <- function(data=NULL,
                          gene.column=NULL,
                          species=NULL){
  if(is.null(data)){
    cat(crayon::red("no data provided, execution stops \n"))
    stop()
  }
  if(class(data)!="data.frame"){
    cat(crayon::red("data should be in data.frame format, execution stops \n"))
    stop()
  }
  if(is.null(gene.column)){
    cat(crayon::red("a column indicating gene IDs must be provided, execution stops \n"))
    stop()
  }
  if(is.null(species)){
    cat(crayon::red("no species provided, execution stops"))
    stop()
  }

  df <- data
  n <- gene.column
  if(species=="human"){k <- matrisome.list$human}
  if(species=="mouse"){k <- matrisome.list$mouse}
  if(species=="c.elegans"){k <- matrisome.list$c.elegans}
  if(species=="zebrafish"){k <- matrisome.list$zebrafish}
  if(species=="drosophila"){k <- matrisome.list$drosophila}
  if(species=="quail"){k <- matrisome.list$quail}

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
    stop()
  }
  if(class(data)!="data.frame"){
    cat(crayon::red("data should be in data.frame format, execution stops \n"))
    stop()
  }
  if(length(attributes(data))<1){
    cat(crayon::red("data should be annotated first, execution stops \n"))
    stop()
  }
  if(attributes(data)$workflow!="matrisomeannotatoR"){
    cat(crayon::red("data should be annotated first, execution stops \n"))
    stop()
  }
  if(ncol(data)<=3){
    cat(crayon::red("input data had only one column and cannot be tabulated, execution stops \n"))
    stop()
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

  d1$Var1 <- factor(d1$Var1,levels = c("Core matrisome",
                                       "Matrisome-associated",
                                       "ECM Glycoproteins",
                                       "Collagens",
                                       "Proteoglycans",
                                       "ECM-affiliated Proteins",
                                       "ECM Regulators",
                                       "Secreted Factors",
                                       "Non-matrisome"
  ))
  nd1 <- data.frame(V1=c(
    "Core matrisome",
    "Matrisome-associated",
    "ECM Glycoproteins",
    "Collagens",
    "Proteoglycans",
    "ECM-affiliated Proteins",
    "ECM Regulators",
    "Secreted Factors",
    "Non-matrisome"),
    V2=c("#002253",
         "#DB3E18",
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

  require(ggalluvial)
  gdf <- data
  d2 <- as.data.frame(table(gdf$`Annotated Matrisome Division`,gdf$`Annotated Matrisome Category`))
  d2 <- d2[d2$Freq>0,]
  d2$Var1 <- factor(d2$Var1,levels=c("Core matrisome",
                                     "Matrisome-associated",
                                     "Non-matrisome"))
  d2$Var2 <- factor(d2$Var2,levels=c("ECM Glycoproteins",
                                     "Collagens",
                                     "Proteoglycans",
                                     "ECM-affiliated Proteins",
                                     "ECM Regulators",
                                     "Secreted Factors",
                                     "Non-matrisome"))

  p2 <- ggplot(d2,aes(y = Freq, axis1 = Var1, axis2 = Var2)) +
    geom_alluvium(aes(fill = Var2), width = 1/12) +
    geom_stratum(width = 1/12, fill = "white", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Annotated Matrisome Division", "Annotated Matrisome Category"), expand = c(.05, .05)) +
    scale_fill_manual(values = c("#13349D",
                                 "#0584B7",
                                 "#59D8E6",
                                 "#F4651E",
                                 "#F9A287",
                                 "#FFE188",
                                 "#D9D9D9")) +
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
#' @param print.plot logical. If TRUE (default) the plot is printed, otherwise the graph object (in its original format) is returned.
#'
#' @return a DonutPie object (printed or not)
#'
#' @examples matripie(data)

matripie <- function(data=NULL,
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
  d2$Var1 <- factor(d2$Var1,levels=c("Core matrisome",
                                     "Matrisome-associated",
                                     "Non-matrisome"))
  d2$Var2 <- factor(d2$Var2,levels=c("ECM Glycoproteins",
                                     "Collagens",
                                     "Proteoglycans",
                                     "ECM-affiliated Proteins",
                                     "ECM Regulators",
                                     "Secreted Factors",
                                     "Non-matrisome"))

  if(isTRUE(print.plot)){
  suppressMessages(
    suppressWarnings(
      PieDonut(d2,aes(Var1,Var2,Freq))
    )
  )}else{
    return(
      PieDonut(d2,aes(Var1,Var2,Freq))
    )
  }
}

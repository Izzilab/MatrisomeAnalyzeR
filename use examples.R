library(MatrisomeAnalyzeR)

#testfile is automatically loaded with package MatrisomeAnalyzeR and comes from the old script in the previous matrisomeDB server

test <- matriannotate(testfile,"EntrezGeneSymbol","human") #example of annotation, just like the app
test2 <- matrianalyze(test) #example of analysis, just like the app

matribar(test) #example of plot, just like the app
# if opting for not printing the plot, the returned object can be further customized:
# p <- matribar(test,print.plot = F)
# library(ggsci)
# p + scale_fill_aaas()

matriflow(test) #extra plotting possibility
matripie(test) #extra plotting possibility


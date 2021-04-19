# Anotación microarray / RNA

## listOfTabs <- DEA_RNASeq
## maPackage <- hgu133a.db
## anotPackage <- org.Hs.eg.db
## ID <- "PROBEID" (para microarray) "ENSEMBL" (para RNASeq)

annotated <- function(listOfTabs, annotPackage, maPackage, ID) {

  for(i in 1:(length(listOfTabs) - 2)) {

    topTab <- listOfTabs[[i]]
    IDtopTab <- cbind(newcol = rownames(topTab), topTab)
    colnames(IDtopTab)[1] <- ID

    if(ID == "PROBEID") {

      myProbes <- rownames(IDtopTab)
      thePackage <- BiocGenerics::eval(parse(text = maPackage))
      geneAnots <- AnnotationDbi::select(thePackage, myProbes,
                                         columns = c("ENTREZID", "SYMBOL", "GENENAME"))

    } else if(ID == "ENSEMBL"){

      geneAnots <- AnnotationDbi::select(annotPackage, keys = rownames(topTab),
                                         keytype = ID,
                                         columns = c("ENSEMBL", "ENTREZID",
                                                     "SYMBOL", "GENENAME"))

    }

    annotTopTab <- merge(x = geneAnots, y = IDtopTab, by.x = ID, by.y = ID)
    annotTopTab <- dplyr::arrange(annotTopTab, adj.P.Val)

    listOfTabs[[i]] <- annotTopTab

  }

  return(listOfTabs)

}


# Analisis de significación biológica microarray / RNA

## listOfTabs
## GOPackage <- org.Hs.egGO
## PATHPackage <- org.Hs.egPATH
## organism <- "human"

bsa <- function(listOfTabs, GOPackage, PATHPackage, organism, pvcoff = 0.05,
                padmethod) {

  listOfSelected <- list()

  for(i in 1:(length(listOfTabs) - 2)) {

    listOfSelected[[i]] <- listOfTabs[[i]]$ENTREZID
    names(listOfSelected)[i] <- names(listOfTables)[i]

  }

  mapped_genes2GO <- AnnotationDbi::mappedkeys(GOPackage)
  mapped_genes2KEGG <- AnnotationDbi::mappedkeys(PATHPackage)
  mapped_genes <- union(mapped_genes2GO , mapped_genes2KEGG)

  comparisonsNames <- names(listOfSelected)

  enrRes <- list()

  for (i in 1:length(listOfSelected)){

    genesIn <- listOfSelected[[i]]
    comparison <- comparisonsNames[i]
    enrich.result <- ReactomePA::enrichPathway(gene = genesIn, organism = organism,
                                               pvalueCutoff = pvcoff,
                                               pAdjustMethod = padmethod,
                                               universe = mapped_genes, readable = TRUE)
    enrRes[[i]] <- enrich.result
    names(enrRes)[i] <- names(listOfSelected)[i]
  }

  return(enrRes)

}

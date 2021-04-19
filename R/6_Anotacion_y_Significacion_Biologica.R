#' Anotación de genes
#'
#' Realiza la anotación de los genes contenidos en un marco de datos, añadiendo
#' a los datos ya almacendados los identificadores de PROBE (para microarrays) o
#' ENSEMBL (para RNA-Seq), de ENTREZ, y los símbolos y nombres de los genes.
#' @param listOfTabs Lista que contiene los datos obtenidos en el análisis de
#'   expresión diferencial. Cada contraste debe estar almacenado en un elemento
#'   de la lista en forma de marco de datos que contenga los nombres de los
#'   probes (microarrays) o de los genes de Ensembl (RNA-Seq) como nombres de
#'   fila. Los dos últimos elementos de la lista no se usan para la anotación.
#' @param maPackage Paquete de anotaciones para microarrays.
#' @param annotPackage Paquete de anotaciones para el organismo estudiado.
#' @param ID Carácter. Identificador de los genes analizados. Para análisis de
#'   microarray \code{ID = "PROBEID"}, mientras que para análisis de RNA-Seq
#'   \code{ID = "ENSEMBL".}
#' @return Transforma la lista introducida añadiendo las anotaciones de los
#'   genes.
#' @export
#' @examples
#' # Para datos de microarrays
#' anot_dea_microarray <- annotated(dea_microarray, maPackage = "hgu133a.db",
#'                                  ID = "PROBEID")
#' head(anot_dea_microarray[[1]])
#'
#' # Para datos de RNA-Seq
#' anot_dea_RNASeq <- annotated(dea_RNASeq, annotPackage = org.Hs.eg.db,
#'                              ID = "ENSEMBL")
#' head(anot_dea_RNASeq[[1]])

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

#' Analisis de significación biológica
#'
#' Realiza el análisis de enriquecimiento mediante el paquete ReactomePA
#' @param listOfTabs Lista resultante de realizar las anotaciones.
#' @param GOPackage Mapas entre los ID de genes Entrez y los ID de Gene Ontology
#'   (GO)
#' @param PATHPackage Mapeos entre identificadores de genes Entrez e
#'   identificadores de vías KEGG
#' @param organism Carácter. Organismo sobre el que se realiza el análisis.
#' @param pvcoff Numérico. Umbral de corte del p-valor para filtrar los
#'   resultados.
#' @param padmethod Carácter. Método utilizado para ajustar los p-valores para
#'   pruebas múltiples. Consulte \code{\link[stats]{p.adjust}} para ver la lista
#'   completa de opciones.
#' @return Una lista con tantos elementos como marcos de datos de genes contenía
#'   la entrada. Cada elemento de la lista es un objeto de la clase
#'   \code{enrichResult}.
#' @export
#' @examples
#' # Para datos de microarrays
#' enRes_microarray <- bsa(anot_dea_microarray, GOPackage = "org.Hs.egGO",
#'                         PATHPackage = "org.Hs.egPATH", organism = "human",
#'                         pvcoff = 0.05, padmethod = "BH")
#' enRes_microarray
#'
#' # Para datos de RNA-Seq
#' enRes_RNASeq <- bsa(anot_dea_RNASeq, GOPackage = "org.Hs.egGO",
#'                     PATHPackage = "org.Hs.egPATH", organism = "human",
#'                     pvcoff = 0.05, padmethod = "BH")
#' enRes_RNASeq

bsa <- function(listOfTabs, GOPackage, PATHPackage, organism, pvcoff = 0.05,
                padmethod) {

  listOfSelected <- list()

  for(i in 1:(length(listOfTabs) - 2)) {

    listOfSelected[[i]] <- listOfTabs[[i]]$ENTREZID
    names(listOfSelected)[i] <- names(listOfTabs)[i]

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

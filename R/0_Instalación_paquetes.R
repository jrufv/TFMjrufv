#' Instalación de paquetes de Bioconductor
#'
#' Instala y carga los paquetes de Bioconductor necesarios para el correcto
#' funcionamiento del paquete en caso de que no lo estuvieran.
#' @param annot_pack Paquetes de anotaciones para los datos del análisis.
#' @return
#' @export
#' @examples
#' ins_pack()

ins_pack <- function(annot_pack = NULL) {
  packBioc <- c("AnnotationDbi", "Biobase", "edgeR", "genefilter", "MSnbase", "oligo",
                "oligoClasses", "pmp", "POMA", "S4Vectors", "SummarizedExperiment",
                "xcms", annot_pack)

  for(i in 1:length(packBioc)) {
    if(!(require(packBioc[i], character.only = TRUE)))
      BiocManager::install(packBioc[i])
    library(packBioc[i], character.only = TRUE)
  }
}

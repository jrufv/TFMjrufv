#' Instalación de paquetes de Bioconductor
#'
#' Instala y carga los paquetes de Bioconductor necesarios para el correcto
#' funcionamiento del paquete en caso de que no lo estuvieran.
#' @param annot_pack Paquetes de anotaciones para los datos del análisis.
#' @return Vacío.
#' @export
#' @examples
#' ins_pack()

ins_pack <- function(annot_pack = NULL) {
  packBioc <- c("AnnotationDbi", "Biobase", "BiocGenerics", "edgeR", "enrichplot",
                "genefilter", "limma", "MSnbase", "oligo", "oligoClasses", "pmp", "POMA",
                "ReactomePA", "S4Vectors", "SummarizedExperiment", "xcms", annot_pack)

  for(i in 1:length(packBioc)) {
    if(!(require(packBioc[i], character.only = TRUE)))
      BiocManager::install(packBioc[i], update = TRUE)
    library(packBioc[i], character.only = TRUE)
  }

  if(!require("TFMjrufv"))
    devtools::install_github("jrufv/TFMjrufv", upgrade = "always")
  library("TFMjrufv")
}

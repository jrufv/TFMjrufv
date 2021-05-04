#' Análisis de Expresión Diferencial
#'
#' Realiza un análisis de expresión diferencial basado en modelos lineales del
#' paquete limma.
#' @param object Objeto de la clase \code{ExpressionSet}
#' @param cont Carácter, especifica los contrastes a realizar. Los nombres de
#'   los grupos experimentales se deben introducir de idéntica manera a como
#'   están codificados en el objeto, separados por un guión.
#' @param name Carácter, especifica nombres alternativos para los contrastes. Si
#'   \code{name = NULL} se usaran nombres iguales a los usados en el parámetro
#'   \code{cont}.
#' @param maxanal Numérico. Número máximo de analitos a mostrar en el resultado.
#' @param adjmethod Carácter. Método utilizado para ajustar los p-valores para pruebas
#'   múltiples. Consulte \code{\link[stats]{p.adjust}} para ver la lista
#'   completa de opciones.
#' @param pvalcoff Numérico. Se filtrarán todos los analitos con un p-valor
#'   mayor al especificado.
#' @param dtmethod Cadena de caracteres que especifica cómo se combinarán los
#'   genes y los contrastes en el esquema de prueba múltiple. Las opciones son
#'   \code{"separate"}, \code{"global"}, \code{"hierarchical"} o
#'   \code{"nestedF"}.
#' @return Una lista con un marco de datos para cada contraste realizado, una
#'   tabla que muestra la cantidad de analitos sobreexpresados e infraexpresados
#'   por contraste y un objeto \code{TestResults} que mantiene los analitos que
#'   han mostrado expresión diferencial en algún contraste. Para cada contraste
#'   se muestra un 1 si está sobreexpresado u -1 si está infraexpresado y un 0
#'   si no hay expresión diferencial. Los marcos de datos de cada contraste
#'   muestran diferente estadísticos obtenidos en el análisis (diferencia media
#'   (logFC), expresión promedio (AveExpr), estadístico t moderado, p-valor,
#'   p-valor ajustado y estadístico B), y se muestra ordenado por p-valor de
#'   forma ascendente.
#' @export
#' @examples
#' # Para los datos de microarrays.
#' dea_microarray <- dea(norm_data_microarray, cont = "Treated-Untreated",
#'                       name = "TRvsUN", maxanal = 1000, adjmethod = "fdr",
#'                       pvalcoff = 0.1, dtmethod = "separate")
#' dea_microarray
#'
#' # Para los datos de RNASeq.
#' dea_RNASeq <- dea(norm_data_RNASeq, cont = "BCell-Kidney", name = "BCvsKi",
#'                   adjmethod = "BH", pvalcoff = 0.1, dtmethod = "global")
#' dea_RNASeq
#'
#' # Para los datos de GC/LC-MS RS.
#' dea_MetabRS <- dea(norm_data_MetabRS, cont = "CD-Control", name = "CDvsCon",
#'                    adjmethod = "holm", pvalcoff = 0.1,
#'                    dtmethod = "hierarchical")
#' dea_MetabRS
#'
#' # Para los datos de contendedores de espectros de MS/NMR (se crea un nuevo
#' # grupo falso para mostrar el funcionamiento con 3 grupos experimentales).
#' newgroups <- c(rep("patient", 15), rep("treated", 15), rep("control", 17))
#' pData(norm_data_MetabSB)[1] <- newgroups
#' dea_MetabSB <- dea(norm_data_MetabSB, cont = c("patient-treated",
#'                    "patient-control", "treated-control"),
#'                    name = c("PAvsTR", "PAvsCO", "TRvsCO"),
#'                    adjmethod = "bonferroni", pvalcoff = 0.1,
#'                    dtmethod = "nestedF")
#' dea_MetabSB
#'
#' # Para los datos de concentraciones de metabolitos (se crean nuevos grupos
#' # falsos para mostrar el funcionamiento de la intersección).
#' newgroups <- c(rep("cac.m", 21), rep("cac.w", 20),
#'                rep("con.m", 15), rep("con.w", 15))
#' pData(norm_data_MetabMC)[1] <- newgroups
#' dea_MetabMC <- dea(norm_data_MetabMC, cont = c("cac.m-con.m", "cac.w-con.w",
#'                    "(cac.m-con.m)-(cac.w-con.w)"), name = c("CACvsCON.M",
#'                    "CACvsCON.W", "INT"), adjmethod = "BY", pvalcoff = 0.1,
#'                    dtmethod = "separate")
#' dea_MetabMC

dea <- function(object, cont = NULL, name = NULL, maxanal = NULL, adjmethod = "BH",
                pvalcoff = NULL, dtmethod = "separate") {

  if(is.null(cont) == TRUE) stop("must define at least one contrast")

  if(is.null(name) == FALSE) {
    if(length(cont) != length(name)){
      stop("you must specify the same number of names as contrasts")
    }
  }

  # Construir matriz de diseño y la matriz de contrastes
  DM <- stats::model.matrix(~ 0 + Group, Biobase::pData(object))
  colnames(DM) <- levels(as.factor(Biobase::pData(object)[[1]]))
  CM <- limma::makeContrasts(contrasts = cont, levels = DM)

  if(is.null(name) == FALSE) {
    colnames(CM) <- name
    cont <- name
  }

  # Identificación de analitos diferencialmente expresados
  if(ncol(Biobase::pData(object)) > 1 & Biobase::varLabels(object)[2] == "lib.size") {

    # Para datos de RNA-Seq
    if(max(Biobase::pData(object)[2]) / min(Biobase::pData(object)[2]) < 3) {
      logCPM <- edgeR::cpm(object, log = TRUE)
      fit <- limma::lmFit(logCPM, DM)
      fit.main <- limma::contrasts.fit(fit, CM)
      istrend = TRUE
    } else {
      v <- limma::voom(object)
      fit <- limma::lmFit(v, DM)
      istrend = FALSE
    }

  } else {

    # Para el resto
    fit <- limma::lmFit(object, DM)
    istrend = FALSE
  }

  fit.main <- limma::contrasts.fit(fit, CM)
  fit.eBayes <- limma::eBayes(fit.main, trend = istrend)

  # Obtener listas de analitos diferencialmente expresados
  if(is.null(maxanal) == TRUE) maxanal <- nrow(fit.eBayes)

  results <- list()

  for(i in 1:length(cont)) {

    topTab <- limma::topTreat(fit.eBayes, coef = cont[i], number = maxanal,
                              adjust.method = adjmethod,
                              p.value = pvalcoff)
    results[[i]] <- topTab

  }

  names(results) <- cont

  # Comparaciones múltiples
  MC <- limma::decideTests(fit.eBayes, method = dtmethod, adjust.method = adjmethod,
                           p.value = pvalcoff)
  sum.MC.rows <- apply(abs(MC), 1, sum)
  MC.selected <- MC[sum.MC.rows != 0, ]

  MC <- summary(MC)
  results[['MultComp']] <- MC
  results[['TestResMat']] <- MC.selected

  return(results)

}

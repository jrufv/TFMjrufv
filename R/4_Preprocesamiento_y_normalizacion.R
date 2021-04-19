#' Preprocesamiento y Normalización de datos de microarrays
#'
#' Realiza el preprocesamiento y la normalización de datos de microarrays
#' almacenados en objetos del tipo \code{ExpressionFeatureSet}.
#' @param object Objeto del tipo \code{ExpressionFeatureSet}.
#' @param varFilter,featureFilter Lógico, si es \code{TRUE} se aplicará el
#'   filtrado por varianza (\code{varFilter}), o el filtrado basado en
#'   anotaciones y la eliminación de duplicados (\code{featureFilter}). Si ambos
#'   son \code{TRUE} se aplican ambos filtrados.
#' @param annot_pack Paquete de anotaciones (para datos de microarrays).
#' @param ... Más argumentos que se pasan a \code{\link[oligo]{rma}} o
#'   \code{\link[genefilter]{nsFilter}}.
#' @return Un objeto de la clase \code{ExpressionSet} que contiene los datos
#'   resultantes del preprocesamiento aplicado
#' @export
#' @examples
#' norm_data_microarray <- prep_norm_microarray(data_microarray,
#'                                              varFilter = TRUE,
#'                                              featureFilter = TRUE,
#'                                              annot_pack = "hgu133a.db")
#' norm_data_microarray

prep_norm_microarray <- function(object, varFilter = TRUE, featureFilter = TRUE,
                                 annot_pack = NULL, ...) {

  dir.create("Results", showWarnings = FALSE)

  if(class(object) != "ExpressionFeatureSet") {
    stop("class object must be ExpressionFeatureSet")
  }

  # Normalización
  norm_obj <- oligo::rma(object)

  # Filtraje no específico
  if(featureFilter == TRUE & base::is.null(annot_pack) == TRUE) {
    stop("must enter an annotations database")
  } else if(featureFilter == TRUE & base::is.null(annot_pack) == FALSE) {
    BiocGenerics::annotation(norm_obj) <- annot_pack
  }

  if(varFilter == TRUE & featureFilter == TRUE) {

    # Filtraje no específico por varianza y anotaciones
    filt_norm_obj <- genefilter::nsFilter(norm_obj)
    data <- filt_norm_obj$eset

  } else if(varFilter == TRUE & featureFilter == FALSE) {

    # Filtraje no específico por varianza
    data <- genefilter::varFilter(norm_obj)

  } else if(varFilter == FALSE & featureFilter == TRUE) {

    # Filtraje no específico por anotaciones
    data <- genefilter::featureFilter(norm_obj)

  } else data <- norm_obj

  # Guardar datos de expresión y objeto
  utils::write.csv(Biobase::exprs(data),
                   file = paste0("./Results/norm_data_microarray_", Sys.time(), ".csv"))

  Biobase::pData(data) <- Biobase::pData(data)[2:ncol(Biobase::pData(data))]

  save(data, file = paste0("./Results/norm_data_microarray_", Sys.time(), ".RData"))

  return(data)

}

#' Preprocesamiento y Normalización de datos de RNA-Seq
#'
#' Realiza el preprocesamiento y la normalización de datos de RNA-Seq
#' almacenados en objetos del tipo \code{DGEList}.
#' @param object Objeto del tipo \code{DGEList}.
#' @param usFilter Lógico, si es \code{TRUE} se aplicará un filtrado no
#'   específico eliminando filas para las que la mitad de las muestras tengan un
#'   recuento por millón menor a 1.
#' @param ... Más argumentos que se pasan a \code{\link[edgeR]{calcNormFactors}}
#' @return Un objeto de la clase \code{DGEList} que contiene los datos
#'   resultantes del preprocesamiento aplicado
#' @export
#' @examples
#' norm_data_RNASeq <- prep_norm_rna(data_RNASeq, usFilter = TRUE,
#'                                   method = "TMM")
#' norm_data_RNASeq

prep_norm_rna <- function(object, usFilter = TRUE, ...) {

  dir.create("Results", showWarnings = FALSE)

  if(class(object) != "DGEList") stop("class object must be DGEList")

  # Filtraje no específico
  if(usFilter == TRUE) {
    keep <- edgeR::filterByExpr(object, min.count = 1,
                                min.total.count = (ncol(object$counts) / 2))
    filt_obj <- object[keep,]
  } else filt_obj <- object

  # Normalización
  data <- edgeR::calcNormFactors(filt_obj)

  # Guardar datos de expresión
  utils::write.csv(data$counts,
                   file = paste0("./Results/norm_data_RNASeq_", Sys.time(), ".csv"))

  # Generar objeto ExpressionSet
  data <- Biobase::ExpressionSet(assayData = data$counts,
                                 phenoData = Biobase::AnnotatedDataFrame(data$samples))
  Biobase::varLabels(data)[1] <- "Group"

  # Guardar objeto
  save(data, file = paste0("./Results/norm_data_RNASeq_", Sys.time(), ".RData"))

  return(data)

}

#' Preprocesamiento de espectros brutos de GC/LC-MS
#'
#' Realiza el preprocesamiento de espectros brutos de GC/LC-MS almacenados en
#' objetos del tipo \code{OnDiskMSnExp} o \code{MSnExp}.
#' @param object Objeto del tipo \code{OnDiskMSnExp} o \code{MSnExp}.
#' @param fCPmethod Método de deteccion de picos a usar. Para más información
#'   visite \code{\link[xcms]{findChromPeaks}}.
#' @param refineRT,refineIn,refineMN Lógico, si es \code{TRUE} se aplica un
#'   filtraje según un determinado tiempo de retención (\code{refineRT}) o
#'   intensidad (\code{refineIn}), o se fusionan picos cromatográficos que se
#'   superponen o presentan unas dimensiones de rt y m/z cercanas
#'   (\code{refineMN}).
#' @param aRtmethod Método de alineamiento de picos a usar. Para más información
#'   visite \code{\link[xcms]{adjustRtime}}.
#' @param gCPmethod Método de correspondencia a usar. Para más información
#'   visite \code{\link[xcms]{groupChromPeaks}}.
#' @param fillCP Lógico, si es \code{TRUE} se realiza un rellenado de los
#'   valores ausentes.
#' @param fillmethod Método de rellenado de valores ausentes a usar. Para más
#'   información visite \code{\link[xcms]{fillChromPeaks}}.
#' @param ... Más argumentos que se pasan a \code{\link[xcms]{findChromPeaks}} o
#'   sus métodos, \code{\link[xcms]{refineChromPeaks}} o sus métodos,
#'   \code{\link[xcms]{adjustRtime}} o sus métodos,
#'   \code{\link[xcms]{groupChromPeaks}} o sus métodos, o
#'   \code{\link[xcms]{fillChromPeaks}} o sus métodos.
#' @return Un objeto de la clase \code{MSnSet} que contiene los datos
#'   resultantes del preprocesamiento aplicado
#' @export
#' @examples
#' prep_data_MetabRS <- prep_metRS(data_MetabRS,
#'                                 fCPmethod = c("centWave", "centWavewpi",
#'                                 "matchedFilter", "massifquant", "MSW"),
#'                                 refineRT = TRUE, refineIn = TRUE,
#'                                 refineMN = TRUE,
#'                                 aRtmethod = c("peakGroups", "obiwarp"),
#'                                 gCPmethod = c("density", "mzClust", "nearest"),
#'                                 fillCP = TRUE, fillmethod = c("fill", "area"))
#' prep_data_MetabRS

prep_metRS <- function(object, fCPmethod = NULL, refineRT = FALSE, refineIn = FALSE,
                       refineMN = FALSE, aRtmethod = NULL, gCPmethod = NULL,
                       fillCP = FALSE, fillmethod = NULL, ...) {

  dir.create("Results", showWarnings = FALSE)

  if(class(object) != "OnDiskMSnExp" & class(object) != "MSnExp") {
    stop("class object must be OnDiskMSnExp or MSnExp")
  }

  if(base::is.null(fCPmethod) == TRUE) {
    stop("must select a chromatographic peak detection method")
  }

  if(base::is.null(aRtmethod) == TRUE) {
    stop("must select an alignment method")
  }

  if(base::is.null(gCPmethod) == TRUE) {
    stop("must select a correspondence method")
  }

  if(fillCP == TRUE & base::is.null(fillmethod) == TRUE) {
    stop("must select a fill method")
  }

  # Detección de picos cromatográficos
  if(fCPmethod == "centWave") {
    find_method <- xcms::CentWaveParam()
  } else if(fCPmethod == "centWavewpi") {
    find_method <- xcms::CentWavePredIsoParam()
  } else if(fCPmethod == "matchedFilter") {
    find_method <- xcms::MatchedFilterParam()
  } else if(fCPmethod == "massifquant") {
    find_method <- xcms::MassifquantParam()
  } else if(fCPmethod == "MSW") {
    find_method <- xcms::MSWParam()
  } else stop("chosen method is not valid")

  find_obj <- xcms::findChromPeaks(object, param = find_method)

  # Refinar la detección de picos
  if(refineRT == TRUE) {
    refine_method <- xcms::CleanPeaksParam()
    rert_find_obj <- xcms::refineChromPeaks(find_obj,
                                            param = refine_method)

  } else rert_find_obj <- find_obj

  if(refineIn == TRUE) {
    refine_method <- xcms::FilterIntensityParam()
    rein_rert_find_obj <- xcms::refineChromPeaks(rert_find_obj,
                                                 param = refine_method)
  } else rein_rert_find_obj <- rert_find_obj

  if(refineMN == TRUE) {
    refine_method <- xcms::MergeNeighboringPeaksParam()
    refi_find_obj <- xcms::refineChromPeaks(rein_rert_find_obj,
                                            param = refine_method)
  } else refi_find_obj <- rein_rert_find_obj

  # Alineamiento de picos
  if(aRtmethod == "peakGroups") {
    adRt_method <- xcms::PeakGroupsParam()
    alig_refi_find_obj <- xcms::adjustRtimePeakGroups(refi_find_obj,
                                                      param = adRt_method)
  } else if(aRtmethod == "obiwarp") {
    adRt_method <- xcms::ObiwarpParam()
    alig_refi_find_obj <- xcms::adjustRtime(refi_find_obj,
                                            param = adRt_method)
  } else stop("chosen method is not valid")

  # Correspondencia
  Groups <- as.character(Biobase::pData(object)[[3]])
  if(gCPmethod == "density") {
    corr_method <- xcms::PeakDensityParam(sampleGroups = Groups)
  } else if(gCPmethod == "nzClust") {
    corr_method <- xcms::MzClustParam(sampleGroups = Groups)
  } else if(gCPmethod == "nearest") {
    corr_method <- xcms::NearestPeaksParam(sampleGroups = Groups)
  } else stop("chosen method is not valid")

  corr_alig_refi_find_obj <- xcms::groupChromPeaks(alig_refi_find_obj,
                                                   param = corr_method)

  # Rellenar picos en valores ausentes
  if(fillCP == TRUE) {
    if(fillmethod == "fill") {
      fill_method <- xcms::FillChromPeaksParam()
    } else if(fillmethod == "area") {
      fill_method <- xcms::ChromPeakAreaParam()
    } else stop("chosen method is not valid")

    data <- xcms::fillChromPeaks(corr_alig_refi_find_obj,
                                 param = fill_method)
  } else data <- corr_alig_refi_find_obj

  # Guardar datos de expresión
  data_save <- xcms::featureValues(data)
  colnames(data_save) <- Biobase::pData(data)[[2]]
  utils::write.csv(data_save,
                   file = paste0("./Results/prep_data_MetabRS_", Sys.time(), ".csv"))

  # Guardar objeto
  save(data, file = paste0("./Results/prep_data_MetabRS_", Sys.time(), ".RData"))

  # Generar objeto MSnSet
  data <- MSnbase::quantify(data)
  data <- MSnbase::as(data, "MSnSet")
  Biobase::pData(data) <- Biobase::pData(data)[3:ncol(Biobase::pData(data))]

  return(data)

}

#' Preprocesamiento de contenedores de espectros de MS o NMR
#'
#' Realiza el preprocesamiento de contenedores de espectros de MS o NMR
#' almacenados en objetos del tipo \code{SummarizedExperiment}.
#' @param object Objeto del tipo \code{SummarizedExperiment}.
#' @param filterMV,filterF,filterRSD,filterB Lógico, si es \code{TRUE} se
#'   realiza el filtraje de los datos. Para \code{filterMV} se filtrarán muestra
#'   en función del porcentaje de valores perdidos, para \code{filterF} se
#'   filtrarán características en función del porcentaje de muestras que
#'   contengan valores no perdidos, para \code{filterRSD} se filtrarán
#'   características en función de la desviación estándar relativa de las
#'   muestras de control de calidad. Para \code{filterB} se filtrarán
#'   características de origen no biológico mediante muestras en blanco.
#' @param mpmv Numérico, valor entre 0 y 1 del umbral del porcentaje de valores
#'   perdidos en la muestra.
#' @param mf Numérico, valor entre 0 y 1 del umbral de fracción de detección.
#' @param fFmethod Método del filtraje de características a usar. Para más
#'   información visite \code{\link[pmp]{filter_peaks_by_fraction}}.
#' @param QClab Carácter, etiqueta de clase utilizada para identificar muestras
#'   de control de calidad. Para \code{filterB = TRUE}, Si el parámetro QClab no
#'   es NULL, se utilizarán muestras de control de calidad para calcular la
#'   intensidad de la señal media.
#' @param mrsd Numérico, valor umbral del %%RSD del control de calidad.
#' @param mfc Numérico, fold change mínimo entre muestras analíticas y en
#'   blanco.
#' @param Blab Carácter, etiqueta de clase utilizada para identificar muestras
#'   en blanco.
#' @param remB Lógico, si es TRUE se eliminarán las muestras en blanco.
#' @param fib Numérico, valor entre 0 y 1 de la fracción de picos en blanco en
#'   las que deben estar presentes.
#' @return Un objeto de la clase \code{MSnSet} que contiene los
#'   datos resultantes del preprocesamiento aplicado
#' @export
#' @examples
#' prep_data_MetabSB <- prep_metSB(data_MetabSB, filterMV = TRUE, mpmv = 0.1,
#'                                 filterF = TRUE, mf = 0.9, fFmethod = "across",
#'                                 filterRSD = FALSE, filterB = FALSE)
#' prep_data_MetabSB

prep_metSB <- function(object, filterMV = TRUE, filterF = TRUE, filterRSD = TRUE,
                       filterB = TRUE, mpmv = NULL, mf = NULL, fFmethod =  NULL,
                       QClab = NULL, mrsd = NULL, mfc = NULL, Blab = NULL,
                       remB = TRUE, fib = 0) {

  dir.create("Results", showWarnings = FALSE)

  if(class(object) != "SummarizedExperiment") {
    stop("class object must be SummarizedExperiment")
  }

  Groups <- as.character(SummarizedExperiment::colData(object)[[1]])

  # Filtrar muestras por valores perdidos
  if(filterMV == TRUE) {

    if(base::is.null(mpmv) == TRUE) {
      stop("must select the maximum allowed percentage of missing values")
    }

    fimv_obj <- pmp::filter_samples_by_mv(object, max_perc_mv = mpmv,
                                          classes = Groups)

  } else fimv_obj <- object

  # Rehacer el objeto Groups por si se han filtrado muestras
  Groups <- as.character(SummarizedExperiment::colData(fimv_obj)[[1]])

  # Filtrar características por valores perdidos
  if(filterF == TRUE) {

    if(base::is.null(mf) == TRUE) {
      stop("must select the minimum allowed relative proportion of samples containing
           non_missing values")
    }

    if(base::is.null(fFmethod) == TRUE) {
      stop("must select a filtered method")
    }

    if(fFmethod != "QC" & fFmethod != "within" & fFmethod != "across") {
      stop("chosen method is not valid")
    }

    if(fFmethod == "QC") {
      if(base::is.null(QClab) == TRUE) {
        stop("must select a quality control class label")
      }
    }

    filF_fimv_obj <- pmp::filter_peaks_by_fraction(fimv_obj, min_frac = mf,
                                                   classes = Groups,
                                                   method = fFmethod,
                                                   qc_label = QClab)

  } else filF_fimv_obj <- fimv_obj

  # Filtrado de características por desviación estándar relativa
  if(filterRSD == TRUE) {

    if(base::is.null(QClab) == TRUE) stop("must select a quality control class label")

    if(is.null(mrsd) == TRUE) {
      stop("must select the maxim allowed threshold of QC RSD% values")
    }

    fRSD_filF_fimv_obj <- pmp::filter_peaks_by_rsd(filF_fimv_obj, max_rsd = mrsd,
                                                   classes = Groups,
                                                   qc_label = QClab)

  } else fRSD_filF_fimv_obj <- filF_fimv_obj

  # Filtrado de características de origen no biológico
  if(filterB == TRUE) {

    if(base::is.null(mfc) == TRUE) {
      stop("must select the minimum fold change between analytical and blanck samples")
    }

    if(is.null(Blab) == TRUE) stop("must select a blank class label")

    data <- pmp::filter_peaks_by_blank(fRSD_filF_fimv_obj, fold_change = mfc,
                                       classes = Groups, blank_label = Blab,
                                       qc_label = QClab, remove_samples = remB,
                                       fraction_in_blank = fib)

  } else data <- fRSD_filF_fimv_obj

  # Guardar datos
  utils::write.csv(SummarizedExperiment::assay(data),
                   file = paste0("./Results/prep_data_MetabSB_", Sys.time(), ".csv"))

  # Guardar objeto
  save(data, file = paste0("./Results/prep_data_MetabSB_", Sys.time(), ".RData"))

  # Generar objeto MSnSet
  data <- MSnbase::as(data, "MSnSet")

  return(data)

}

#' Imputación y normalización de datos metabolómicos
#'
#' Realiza la imputación de valores perdidos, normalización y detección y
#' eliminación de outliers a datos metabolómicos en objetos del tipo
#' \code{MSnSet}.
#' @param object Objeto del tipo \code{MSnSet}.
#' @param impute,routliers Lógico, si es \code{TRUE} se realizará la imputación
#'   de valores perdidos (para \code{impute} o la detección y eliminación de
#'   outliers (para \code{routliers}).
#' @param coff Numérico, indica el porcentaje de valores perdidos permitido en
#'   cada grupo. Si uno de los grupos tiene menos valores perdidos que el valor
#'   de corte seleccionado, esta característica no se eliminará.
#' @param immethod Método de imputación de valores perdidos a usar. Las opciones
#'   son: \code{"none"}, \code{"half_min"}, \code{"median"}, \code{"mean"},
#'   \code{"min"}, \code{"knn"} y \code{"rf"}. Si es \code{"none"}
#'   (predeterminado), todos los valores perdidos serán reemplazados por cero.
#' @param nomethod Método de normalización. Las opciones son: \code{"none"}
#'   (predeterminado), \code{"auto_scaling"}, \code{"level_scaling"},
#'   \code{"log_scaling"}, \code{"log_transformation"}, \code{"vast_scaling"} y
#'   \code{"log_pareto"}.
#' @param oumethod Método de detección de outliers. Las opciones son
#'   \code{"median"} y \code{"centroid"}.
#' @param ditype Tipo de medida de distancia. Las opciones son
#'   \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}
#'   y \code{"minkiowski"}.Vea \code{\link[stats]{dist}}.
#' @param ... Más argumentos que se pasan a \code{\link[POMA]{PomaImpute}}.
#' @return Un objeto de la clase \code{MSnSet} que contiene los datos
#'   resultantes del preprocesamiento aplicado
#' @export
#' @examples
#' norm_data_MetabMC <- met_imp_norm(data_MetabMC, impute = TRUE, coff = 20,
#'                                   immethod = "knn", nomethod = "log_pareto",
#'                                   routliers = TRUE, oumethod = "median",
#'                                   ditype = "euclidean")
#' norm_data_MetabMC

met_imp_norm <- function(object, impute = TRUE, coff = 20, immethod = "none",
                         nomethod = "none", routliers = TRUE, oumethod = NULL,
                         ditype = "euclidean", ...) {

  dir.create("Results", showWarnings = FALSE)

  if(class(object) != "MSnSet") stop("class object must be MSnSet")

  # Imputación de valores perdidos
  if(impute == TRUE) {

    impu_obj <- POMA::PomaImpute(object, cutoff = coff, method = immethod)

  } else impu_obj <- object

  # Normalización
  norm_impu_obj <- POMA::PomaNorm(impu_obj, method = nomethod)

  # Detección y eliminación de outliers
  if(routliers == TRUE) {

    if(base::is.null(oumethod) == TRUE) stop("must select a remove outliers method")

    if(oumethod != "median" & oumethod != "centroid") {
      stop("chosen method is not valid")
    }

    if(ditype != "euclidean" & ditype != "maximum" & ditype != "manhattan" &
       ditype != "canberra" & ditype != "minkowski") {
      stop("chosen distance measure method is not valid")
    }

    data <- POMA::PomaOutliers(norm_impu_obj, do = "clean", method = ditype,
                               type = oumethod)

  } else data <- norm_impu_obj


  # Guardar datos
  utils::write.csv(Biobase::exprs(data),
                   file = paste0("./Results/norm_Metab_", Sys.time(), ".csv"))

  # Guardar objeto
  save(data, file = paste0("./Results/norm_Metab_", Sys.time(), ".RData"))

  # Generar objeto ExpressionSet
  data <- MSnbase::as(data, "ExpressionSet")

  return(data)

}

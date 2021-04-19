#' Lectura de archivos de datos ómicos
#'
#' Lectura de diferentes tipos de archivos de datos ómicos y construcción de un
#' objeto en el que almacenar los datos de los ensayos y las muestras.
#' @param data_type Tipo de datos a introducir, puede ser \code{"microarray"},
#'   \code{"RNA-Seq"}, \code{"MetabRS"}, \code{"MetabSB"} o \code{"MetabMC"}.
#' @param path Ruta de ubicación de los archivos para \code{data_type =
#'   "microarray"} o \code{"MetabRS"}.
#' @param file_type Tipo de archivos, para \code{data_type = "MetabRS"}. Puede
#'   ser \code{".NetCDF"}, \code{".mzML"}, \code{".mzXML"} o \code{".mzData"}.
#' @param raw_data Ruta y nombre del archivo con los datos brutos del
#'   experimento, para \code{data_type = "RNA-Seq"}, \code{"MetabSB"} o
#'   {"MetabMC"}. Debe ser un archivo \code{.csv} o \code{.txt} con cabecera. La
#'   primera columna debe corresponder a los genes (\code{"RNA-Seq"}), los bins
#'   (\code{"MetabSB"}) o los metabolitos (\code{"MetabMC"}) y las siguientes a
#'   las muestras.
#' @param sep_rd Carácter separador de campo para el archivo \code{raw_data}.
#'   Los valores de cada línea del archivo están separados por este carácter. De
#'   forma predeterminada \code{sep_rd = ""}.
#' @param targets Ruta y nombre del archivo con la información sobre las
#'   muestras del experimento. Debe ser un archivo \code{.csv} o \code{.txt} con
#'   cabecera. Para \code{data_type = "microarray"} o \code{"MetabRS"} debe
#'   tener tres columnas como mínimo, la primera con el nombre de los archivos
#'   que contienen los datos de cada muestra, la segunda con un identificador
#'   único de las muestras y la tercera con el grupo experimental de las
#'   muestras. Para el resto debe tener dos columnas como mínimo, la primera con
#'   el identificador de las muestras y la segunda con el grupo experimental de
#'   las muestras. En ambos casos se aceptan más columnas como covariables.
#' @param sep_targ Carácter separador de campo para el archivo \code{targets}.
#'   Los valores de cada línea del archivo están separados por este carácter. De
#'   forma predeterminada \code{sep_targ = ""}.
#' @param mode Puede ser \code{"inMemory"} (predeterminado) o \code{"onDisk"}.
#'   El primero carga los datos sin procesar en la memoria, mientras que el
#'   segundo solo genera el objeto y se accede a los datos sin procesar en el
#'   disco cuando es necesario (sólo válido para \code{data_type = "MetabRS"}).
#' @return Un objeto de la clase \code{ExpressionFeatureSet} (para microarrays),
#'   \code{DGEList} (para RNA-Seq), \code{OnDiskMSnExp} o \code{MSnExp} (para
#'   espectros brutos de GC/LC-MS), \code{SummarizedExperiment} (para
#'   contenedores de espectros de MS/NMR), o \code{MSnSet} (para concentraciones
#'   de metabolitos) que contiene los datos de los ensayos y datos sobre las
#'   muestras.
#' @export
#' @examples
#' # Para microarrays
#' data_microarray <- read_data(data_type = "microarray",
#'                              path = "./data/microarray",
#'                              targets = "./data/microarray/targets.csv",
#'                              sep_targ = ";")
#' data_microarray
#'
#' # Para RNA-Seq
#' data_RNASeq <- read_data(data_type = "RNA-Seq",
#'                          raw_data = "./data/RNA-Seq/counts.csv",
#'                          sep_rd = "\t",
#'                          targets = "./data/RNA-Seq/targets.csv",
#'                          sep_targ = ",")
#' data_RNASeq
#'
#' # Para Espectros Brutos de GC/LC-MS
#' data_MetabRS <- read_data(data_type = "MetabRS",
#'                           path = "./data/Met_mzML",
#'                           file_type = ".mzML",
#'                           targets = "./data/Met_mzML/targets.csv",
#'                           sep_targ = ";",
#'                           mode = "onDisk")
#' data_MetabRS
#'
#' # Para Contenedores de Espectros de MS/NMR
#' data_MetabSB <- read_data(data_type = "MetabSB",
#'                           raw_data = "./data/Met_spectra_bins/spectra.csv",
#'                           sep_rd = ",",
#'                           targets = "./data/Met_spectra_bins/targets.csv",
#'                           sep_targ = ",")
#' data_MetabSB
#'
#' # Para Concentraciones de Metabolitos
#' data_MetabMC <- read_data(data_type = "MetabMC",
#'                          raw_data = "./data/Met_conc/concent.csv",
#'                          sep_rd = ",",
#'                          targets = "./data/Met_conc/targets.csv",
#'                          sep_targ = ",")
#' data_MetabMC

read_data <- function(data_type, path, file_type, raw_data, sep_rd = "",
                      targets, sep_targ = "", mode = "onDisk") {

  if(missing(data_type)) stop("argument data_type is missing, with no default")
  if(data_type != "microarray" & data_type != "RNA-Seq" & data_type != "MetabRS" &
     data_type != "MetabPL" & data_type != "MetabSB" & data_type != "MetabMC") {
    stop("Data type must be microarray, RNA-Seq, MetabRS, MetabPL, MetabSB or MetabMC")
  }

  if(data_type == "microarray" | data_type == "MetabRS" | data_type == "MetabPL") {
    if(missing(path)) stop("argument path is missing, with no default")

    min_col = 3
  }

  if(data_type == "MetabRS") {
    if(missing(file_type)) stop("argument file_type is missing, with no default")
    if(file_type != ".NetCDF" & file_type != ".mzML" &
       file_type != ".mzXML" & file_type != ".mzData") stop("the file_type is wrong")
  }

  if(data_type == "RNA-Seq" | data_type == "MetabSB" | data_type == "MetabMC") {
    if(missing(raw_data)) stop("argument counts is missing, with no default")

    file <- utils::read.table(raw_data, header = TRUE, row.names = 1, sep = sep_rd)

    min_col = 2
  }

  if(missing(targets)) stop("argument targets is missing, with no default")

  if(data_type != "microarray") {
    sampleinfo <- utils::read.table(targets, header = TRUE, sep = sep_targ,
                                    stringsAsFactors = TRUE)
    if(ncol(sampleinfo) < min_col) {
      stop("Insufficient sample information")
    }
  }

  if(data_type == "microarray") {

    files <- oligoClasses::list.celfiles(path, full.names = TRUE)
    sampleinfo <- Biobase::read.AnnotatedDataFrame(targets, header = TRUE,
                                                   row.names = 1, sep = sep_targ)
    min_col <- 2
    if(ncol(sampleinfo) < min_col) {
      stop("Insufficient sample information")
    }
    data <- oligo::read.celfiles(files, phenoData = sampleinfo)
    sampleinfo@data[,1] -> rownames(Biobase::pData(data))
    colnames(data) <- rownames(Biobase::pData(data))

  } else if(data_type == "RNA-Seq") {

    if(ncol(sampleinfo) == min_col) {
      data <- edgeR::DGEList(file, group = sampleinfo[,2])
    } else {
      data <- edgeR::DGEList(file, group = sampleinfo[,2],
                             samples = sampleinfo[,3:ncol(sampleinfo)])
    }

  } else if(data_type == "MetabRS") {

    files <- list.files(path, recursive = TRUE, full.names = TRUE, pattern = file_type)
    data <- MSnbase::readMSData(files = files,
                                pdata = methods::new("NAnnotatedDataFrame", sampleinfo),
                                mode = mode)
    sampleinfo[,2] -> rownames(Biobase::pData(data))

  } else if(data_type == "MetabPL") {

    stop("en contstrucción")

  } else if(data_type == "MetabSB") {

    rdata <- S4Vectors::DataFrame(row.names = rownames(file))
    cdata <- S4Vectors::DataFrame(sampleinfo[2:ncol(sampleinfo)],
                                  row.names = sampleinfo[,1])
    file <- list(file)
    file <- S4Vectors::SimpleList(file)
    data <- SummarizedExperiment::SummarizedExperiment(file, rowData = rdata,
                                                       colData = cdata)

  } else if(data_type == "MetabMC") {

    file <- t(file)
    data <- POMA::PomaMSnSetClass(target = sampleinfo, features = file)

  }

  return(data)

}

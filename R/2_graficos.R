#' Gráficos
#'
#' Realiza diferentes tipos de gráficos en función del objeto introducido.
#' Acepta objetos del tipo \code{ExpressionFeatureSet}, \code{ExpressionSet},
#' \code{DGEList}, \code{OnDiskMSnExp}, \code{MSnExp},
#' \code{SummarizedExperiment} o \code{MSnSet}.
#' @param object Objeto del tipo \code{ExpressionFeatureSet},
#'   \code{ExpressionSet}, \code{DGEList}, \code{OnDiskMSnExp}, \code{MSnExp},
#'   \code{SummarizedExperiment} o \code{MSnSet}.
#' @param plot Tipo de gráfico. Puede ser \code{"PCA"} para gráfico de
#'   componentes principales, \code{"boxplot"} para diagramas de caja,
#'   \code{"heatmap"} para un mapa de calor de la varianza de genes,
#'   \code{"barplot"} para un gráfico de barras de las librerías usadas (sólo
#'   para datos de RNA-Seq), \code{"BPC"} para un cromatograma de pico base
#'   (sólo para espectros de GC/LC-MS), \code{"heatmapS"} para un mapa de calor
#'   por muestras o \code{"density"} para un gráfico de densidad por muestra.
#' @return Un gráfico distinto en función de la selección realizada.
#' @export
#' @examples
#' plotTFM(object, plot = c("PCA", "boxplot", "heatmap", "barplot", "BPC",
#'                          "heatmapS", "density")

plotTFM <- function (object, plot) {

  if(class(object) == "ExpressionFeatureSet") {
    dat_mat <- Biobase::exprs(object)
    Groups <- Biobase::pData(object)[[2]]
    labels <- Biobase::sampleNames(object)
    name <- Biobase::varLabels(object)[2]
    lev <- levels(as.factor(Biobase::pData(object)[[2]]))
    title <- "Microarrays"
    value <- "Expression"
  } else if(class(object) == "ExpressionSet") {
    dat_mat <- Biobase::exprs(object)
    Groups <- Biobase::pData(object)[[1]]
    labels <- Biobase::sampleNames(object)
    name <- Biobase::varLabels(object)[1]
    lev <- levels(as.factor(Biobase::pData(object)[[1]]))
    title <- "Microarrays"
    value <- "Expression"
  } else if(class(object) == "DGEList") {
    dat_mat <- object$counts
    Groups <- as.character(object$samples[[1]])
    labels <- colnames(object$counts)
    name <- colnames(object$samples)[1]
    lev <- levels(object$samples[[1]])
    title <- "RNA-Seq"
    value <- "Counts"
  } else if(class(object) == "OnDiskMSnExp" | class(object) ==  "MSnExp") {
    dat_mat <- matrix(MSnbase::tic(object), ncol = nrow(Biobase::pData(object)))
    colnames(dat_mat) <- rownames(Biobase::pData(object))
    Groups <- as.character(Biobase::pData(object)[[3]])
    rownames(dat_mat) <- Biobase::featureNames(
      Biobase::featureData(object))[1:(length(MSnbase::tic(object))/length(Groups))]
    labels <- rownames(Biobase::pData(object))
    name <- as.character(Biobase::pData(object)[[2]])
    lev <- levels(Biobase::pData(object)[[3]])
    title <- "Raw Spectra"
    value <- "Intensity"
  } else if(class(object) == "SummarizedExperiment") {
    dat_mat <- as.matrix(SummarizedExperiment::assay(object))
    Groups <- as.character(SummarizedExperiment::colData(object)[[1]])
    labels <- colnames(object)
    name <- colnames(SummarizedExperiment::colData(object))
    lev <- levels(SummarizedExperiment::colData(object)[[1]])
    title <- "Spectra Bins"
    value <- "Intensity"
  } else if(class(object) == "MSnSet") {
    dat_mat <- Biobase::exprs(object)
    Groups <- as.character(Biobase::pData(object)[[1]])
    labels <- Biobase::sampleNames(object)
    name <- Biobase::varLabels(object)
    lev <- levels(Biobase::pData(object)[[1]])
    title <- "Metabolite Concentration"
    value <- "Concentration"
  }

  if(class(object) == "ExpressionSet") {
    if(min(Biobase::exprs(object)) < 0) {
      datacpm <- dat_mat
    } else {
      datacpm <- edgeR::cpm(dat_mat, log = TRUE)
    }
  } else {
    datacpm <- edgeR::cpm(dat_mat, log = TRUE)
  }

  meltdatacpm <- reshape::melt(datacpm)
  meltdatacpm <- data.frame(meltdatacpm, condition = cond_vect(object))
  meltdatacpm <- meltdatacpm[,-1]
  names(meltdatacpm) <- c("samples", "value", "condition")

  meltdata <- reshape::melt(dat_mat)
  meltdata <- data.frame(meltdata, condition = cond_vect(object))
  meltdata <- meltdata[,-1]
  names(meltdata) <- c("samples", "value", "condition")

  if(plot == "PCA") {

    pcmat <- stats::prcomp(t(dat_mat), scale=FALSE)
    pcdata <- data.frame(pcmat$x)
    loads <- round(pcmat$sdev^2/sum(pcmat$sdev^2)*100,1)

    ggplot2::ggplot(pcdata, ggplot2::aes(x = pcdata$PC1, y = pcdata$PC2)) +
      ggplot2::theme_bw() +
      ggplot2::geom_hline(yintercept = 0, color = "gray70") +
      ggplot2::geom_vline(xintercept = 0, color = "gray70") +
      ggplot2::geom_point(ggplot2::aes(color = Groups)) +
      ggplot2::coord_cartesian(xlim = c(min(pcmat$x[,1]) - 5, max(pcmat$x[, 1]) + 5)) +
      ggrepel::geom_text_repel(ggplot2::aes(label = labels), size = 3) +
      ggplot2::labs(x = c(paste("PC1", loads[1], "%")),
                    y = c(paste("PC2", loads[2], "%"))) +
      ggplot2::ggtitle(paste("PCA para: ", title, sep=" ")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::scale_color_manual(values = grDevices::topo.colors(ngroup(object)))

  } else if(plot == "boxplot") {

    ggplot2::ggplot(meltdatacpm, ggplot2::aes(x = meltdatacpm$samples, y = value,
                                              fill = meltdatacpm$condition)) +
      ggplot2::theme_bw() +
      ggplot2::geom_boxplot(notch = TRUE) +
      ggplot2::labs(x = "Samples", y = paste(value, "(cpm)", sep = " ")) +
      ggplot2::ggtitle(paste("Boxplot para:", title, sep = " ")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::scale_fill_manual(values = grDevices::topo.colors(ngroup(object))) +
      ggplot2::scale_x_discrete(limits = labels)

  } else if(plot == "heatmap") {

    var_rows <- apply(datacpm, 1, stats::var)
    if(nrow(datacpm) > 1000) n = 500 else n = nrow(datacpm) * 0.25
    select_var <- names(sort(var_rows, decreasing = TRUE))[1:n]
    datacpm <- datacpm[select_var,]
    mypalette <- RColorBrewer::brewer.pal(11,"RdYlBu")
    morecols <- grDevices::colorRampPalette(mypalette)

    gplots::heatmap.2(datacpm,  scale = "row", col = rev(morecols(50)), trace = "none",
                      ColSideColors = col_vect(object), cexCol = 1, srtCol = 0,
                      adjCol = 0.5, main = paste("Heatmap para:", title, sep = " "),
                      xlab = "Samples", ylab = "Compounds")

  } else if(plot == "barplot") {

    if(class(object) == "DGEList") {

      datalibsize <- data.frame(samples = labels, lib.size = object$samples$lib.size,
                                group = Groups)

      ggplot2::ggplot(datalibsize,
                      ggplot2::aes(x = datalibsize$samples,y = datalibsize$lib.size,
                                   fill = datalibsize$group)) +
        ggplot2::theme_bw() +
        ggplot2::geom_bar(stat = "identity", width = 0.5) +
        ggplot2::geom_text(ggplot2::aes(label = datalibsize$lib.size),vjust = 1.5,
                           color = "white", size = 3.5) +
        ggplot2::labs(x = "Samples", y = "Library size") +
        ggplot2::ggtitle("Tamaño de librería por muestra para: RNA-Seq") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::scale_fill_manual(values = grDevices::topo.colors(ngroup(object)))

    } else stop("Barplot is only valid for DGEList objects")

  } else if(plot == "BPC") {

    if(class(object) == "OnDiskMSnExp" | class(object) ==  "MSnExp") {

      bpis <- MSnbase::chromatogram(object, aggregationFun = "max")
      MSnbase::plot(bpis, col = col_vect(object))
      graphics::legend("topright", fill = grDevices::topo.colors(ngroup(object)),
                       legend = lev)

    } else stop("BPC is only valid for OnDiskMSnExp objects")

  } else if(plot == "heatmapS") {

    cormat <- stats::cor(dat_mat)
    mypalette <- RColorBrewer::brewer.pal(11,"RdYlBu")
    morecols <- grDevices::colorRampPalette(mypalette)

    gplots::heatmap.2(cormat, col = rev(morecols(50)), trace = "none",
                      ColSideColors = col_vect(object), cexCol = 1, cexRow = 1,
                      srtCol = 0, adjCol = 0.5,
                      main = paste("Heatmap de muestras para:", title, sep = " "),
                      xlab = "Samples", dendrogram = "col")
    graphics::legend("bottomleft", fill = grDevices::topo.colors(ngroup(object)), legend = lev)

  } else if(plot == "density") {

    ggplot2::ggplot(meltdata, ggplot2::aes(x = value, color = meltdata$condition)) +
      ggplot2::theme_bw() +
      ggplot2::geom_density() +
      ggplot2::labs(x = value, y = "Density") +
      ggplot2::ggtitle(paste("Gráfico de densidad para:", title, sep = " ")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::scale_color_manual(values = grDevices::topo.colors(ngroup(object)))

  }

}

#' Vector de colores
#'
#' Genera un vector de colores por muestra para objetos
#' \code{ExpressionFeatureSet}, \code{ExpressionSet}, \code{DGEList},
#' \code{OnDiskMSnExp}, \code{MSnExp}, \code{SummarizedExperiment} y
#' \code{MSnSet}.
#' @param object Objeto.
#' @return vector de caracteres (colores).
#' @export
#' @examples
#' col_vect(data_microarray)

col_vect <- function(object) {
  if(class(object) == "ExpressionFeatureSet") {
    vector <- Biobase::pData(object)[[2]]
    factor <- as.factor(vector)
  } else if(class(object) == "ExpressionSet") {
    vector <- as.character(Biobase::pData(object)[[1]])
    factor <- as.factor(Biobase::pData(object)[[1]])
  } else if(class(object) == "DGEList") {
    factor <- object$samples$group
    vector <- as.character(factor)
  } else if(class(object) == "OnDiskMSnExp" | class(object) ==  "MSnExp") {
    factor <- Biobase::pData(object)[[3]]
    vector <- as.character(factor)
  } else if(class(object) == "SummarizedExperiment") {
    factor <- SummarizedExperiment::colData(object)[[1]]
    vector <- as.character(factor)
  } else if(class(object) == "MSnSet") {
    factor <- Biobase::pData(object)[[1]]
    vector <- as.character(factor)
  }

  colors <- grDevices::topo.colors(length(levels(factor)))

  for(i in 1:length(levels(factor))) {
    for(j in 1:length(vector)) {
      if(levels(factor)[i] == vector[j]) vector[j] <- colors[i]
    }
  }
  return(vector)
}

#' Vector de condición
#'
#' Crear un vector con el grupo experimental para cada muestra y compuesto a
#' partir de objetos \code{ExpressionFeatureSet}, \code{ExpressionSet},
#' \code{DGEList}, \code{OnDiskMSnExp}, \code{MSnExp},
#' \code{SummarizedExperiment} o \code{MSnSet}.
#' @param object Objeto.
#' @return Vector de caracteres (condición experimental).
#' @export
#' @examples
#' cond_vect(data_microarray)

cond_vect <- function(object) {

  if(class(object) == "ExpressionFeatureSet") {
    datos <- Biobase::exprs(object)
    grupos <- Biobase::pData(object)[[2]]
  } else if(class(object) == "ExpressionSet") {
    datos <- Biobase::exprs(object)
    grupos <- as.character(Biobase::pData(object)[[1]])
  } else if(class(object) == "DGEList") {
    datos <- object$counts
    grupos <- as.character(object$samples$group)
  } else if(class(object) == "OnDiskMSnExp" | class(object) ==  "MSnExp") {
    datos <- matrix(MSnbase::tic(object), ncol = nrow(Biobase::pData(object)))
    colnames(datos) <- rownames(Biobase::pData(object))
    grupos <- as.character(Biobase::pData(object)[[3]])
  } else if(class(object) == "SummarizedExperiment") {
    datos <- SummarizedExperiment::assay(object)
    grupos <- as.character(SummarizedExperiment::colData(object)[[1]])
  } else if(class(object) == "MSnSet") {
    datos <- Biobase::exprs(object)
    grupos <- as.character(Biobase::pData(object)[[1]])
  }

  condition <- c()
  for(i in 1:ncol(datos)) {
    condition <- c(condition, rep(grupos[i], nrow(datos)))
  }
  return(condition)
}

#' Grupos experimentales
#'
#' Devuelve el número de grupos experimentales de objetos
#' \code{ExpressionFeatureSet}, \code{ExpressionSet}, \code{DGEList},
#' \code{OnDiskMSnExp}, \code{MSnExp}, \code{SummarizedExperiment} y
#' \code{MSnSet}.
#' @param object Objeto.
#' @return Numérico (número de grupos experimentales).
#' @export
#' @examples
#' ngroup(data_microarray)

ngroup <- function(object) {
  if(class(object) == "ExpressionFeatureSet") {
    groups <- length(levels(as.factor(Biobase::pData(object)[[2]])))
  } else if(class(object) == "ExpressionSet") {
    groups <- length(levels(as.factor(Biobase::pData(object)[[1]])))
  } else if(class(object) == "DGEList") {
    groups <- length(levels(object$samples[[1]]))
  } else if(class(object) == "OnDiskMSnExp" | class(object) ==  "MSnExp") {
    groups <- length(levels(Biobase::pData(object)[[3]]))
  } else if(class(object) == "SummarizedExperiment") {
    groups <- length(levels(as.data.frame(SummarizedExperiment::colData(object))[[1]]))
  } else if(class(object) == "MSnSet") {
    groups <- length(levels(Biobase::pData(object)[[1]]))
  }
  if(groups == 0) groups == 1
  return(groups)
}

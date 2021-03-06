% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_graficos.R
\name{cond_vect}
\alias{cond_vect}
\title{Vector de condición}
\usage{
cond_vect(object)
}
\arguments{
\item{object}{Objeto.}
}
\value{
Vector de caracteres (condición experimental).
}
\description{
Crear un vector con el grupo experimental para cada muestra y compuesto a
partir de objetos \code{ExpressionFeatureSet}, \code{ExpressionSet},
\code{DGEList}, \code{OnDiskMSnExp}, \code{MSnExp},
\code{SummarizedExperiment} o \code{MSnSet}.
}
\examples{
cond_vect(data_microarray)
}

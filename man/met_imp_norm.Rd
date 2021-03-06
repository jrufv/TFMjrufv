% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_Preprocesamiento_y_normalizacion.R
\name{met_imp_norm}
\alias{met_imp_norm}
\title{Imputación y normalización de datos metabolómicos}
\usage{
met_imp_norm(
  object,
  impute = TRUE,
  coff = 20,
  immethod = "none",
  nomethod = "none",
  routliers = TRUE,
  oumethod = NULL,
  ditype = "euclidean",
  ...
)
}
\arguments{
\item{object}{Objeto del tipo \code{MSnSet}.}

\item{impute, routliers}{Lógico, si es \code{TRUE} se realizará la imputación
de valores perdidos (para \code{impute} o la detección y eliminación de
outliers (para \code{routliers}).}

\item{coff}{Numérico, indica el porcentaje de valores perdidos permitido en
cada grupo. Si uno de los grupos tiene menos valores perdidos que el valor
de corte seleccionado, esta característica no se eliminará.}

\item{immethod}{Método de imputación de valores perdidos a usar. Las opciones
son: \code{"none"}, \code{"half_min"}, \code{"median"}, \code{"mean"},
\code{"min"}, \code{"knn"} y \code{"rf"}. Si es \code{"none"}
(predeterminado), todos los valores perdidos serán reemplazados por cero.}

\item{nomethod}{Método de normalización. Las opciones son: \code{"none"}
(predeterminado), \code{"auto_scaling"}, \code{"level_scaling"},
\code{"log_scaling"}, \code{"log_transformation"}, \code{"vast_scaling"} y
\code{"log_pareto"}.}

\item{oumethod}{Método de detección de outliers. Las opciones son
\code{"median"} y \code{"centroid"}.}

\item{ditype}{Tipo de medida de distancia. Las opciones son
\code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}
y \code{"minkiowski"}.Vea \code{\link[stats]{dist}}.}

\item{...}{Más argumentos que se pasan a \code{\link[POMA]{PomaImpute}}.}
}
\value{
Un objeto de la clase \code{MSnSet} que contiene los datos
resultantes del preprocesamiento aplicado
}
\description{
Realiza la imputación de valores perdidos, normalización y detección y
eliminación de outliers a datos metabolómicos en objetos del tipo
\code{MSnSet}.
}
\examples{
norm_data_MetabMC <- met_imp_norm(data_MetabMC, impute = TRUE, coff = 20,
                                  immethod = "knn", nomethod = "log_pareto",
                                  routliers = TRUE, oumethod = "median",
                                  ditype = "euclidean")
norm_data_MetabMC
}

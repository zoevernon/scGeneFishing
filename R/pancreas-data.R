#' scRNA seq matrix of pancreas cells
#'
#' Data from Tabula Muris project.  A subset of the single cell RNA-seq data 
#' from the pancreas to be used as an example for GeneFishing.  The rows are 
#' genes and columns are cells.
#'
#' @docType data
#'
#' @usage data(pancreas)
#'
#' @keywords datasets
#'
#' @references Schaum, N. et al. (2018) Nature 562, 367–372
#' (\href{https://www.nature.com/articles/s41586-018-0590-4}{Nature})
#'
#' @source \href{https://tabula-muris.ds.czbiohub.org/}{Tabula Muris database}
#'
#' @examples
#' data(pancreas)
#' pancreas[1:4, 1:4]
"pancreas"
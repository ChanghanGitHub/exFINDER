#' calculate average expression levels of input genes in every cluster
#'
#' @param Gene user-select genes (a list of gene symbols)
#' @param Exp.Data expression matrix (rows: features, columns: cells)
#' @param Meta.Data meta data of the expression matrix, must have a variable "Type" containing all cell labels
#'
#' @return
#' @export
#'
#' @examples
get_AveExp <- function(Gene,
                       Exp.Data,
                       Meta.Data){

  Cellgroup <- levels(Meta.Data$Type)
  m <- length(Cellgroup)

  output = data.frame(Gene)

  for (i in 1:m) {

    cell.use <- rownames(Meta.Data)[Meta.Data$Type %in% Cellgroup[i]]

    if (length(cell.use)==0){
      cat("No cell in cell type", Cellgroup[i] ,"\n")
      x <- rep("",each = length(Gene))
      output[paste(Cellgroup[i])] <- x
    } else {
      loc <- match(Gene, row.names(Exp.Data))
      Exp <- as.data.frame(Exp.Data[loc, cell.use])
      if (length(Gene)==1){
        Exp <- as.data.frame(t(Exp))
      }
      n <- length(Exp[1, ])
      Ave.Exp <- rowSums(Exp[, 1:n])/n
      output[paste(Cellgroup[i])] <- Ave.Exp
    }

  }

  return(output)
}

#' select the highest expression level among all clusters based on the average expression levels (data.frame)
#'
#' @param AveExp dataframe of averagen expression levels of genes in clusters, output of "get_AveExp"
#' @importFrom stringr str_c
#'
#' @return
#' @export
#'
#' @examples
get_RowMax <- function(AveExp){

  n <- nrow(AveExp)
  m <- ncol(AveExp)
  output.2 <- data.frame()
  Gene <- AveExp[, 1]



  if (m==2){
    RowData <- AveExp
    out <- names(AveExp)[2] %>%
      str_c(collapse = ",")
    RowData$Max_Cellgroup <- out
    RowData$Max_AveExp <- RowData[, 2]
    output.2 <- rbind(output.2, RowData)
  } else if (m>2){
    AveExp <- AveExp[, -1]
    for (i in 1:n) {
      RowData <- AveExp[i, ]
      index <- which(RowData == max(RowData))
      if (length(index)>1){
        index <- index[1]
      }
      out <- names(RowData[index]) %>%
        str_c(collapse = ",")
      max <- as.numeric(RowData[index])

      RowData$Max_Cellgroup <- out
      RowData$Max_AveExp <- max
      output.2 <- rbind(output.2, RowData)
    }
    output.2 <- cbind(Gene, output.2)
  }

  return(output.2)
}

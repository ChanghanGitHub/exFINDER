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
#' @importFrom magrittr %>%
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

#' calculate the precentiel values of the cell groups in a dataset
#'
#' @param Exp.Data expression matrix (rows: features, columns: cells)
#' @param Meta.Data meta data of the expression matrix, must have a variable "Type" containing all cell labels
#' @param percentile a vector of probabilities (default is c(.5, .75, .90))
#'
#' @return
#' @export
#'
#' @examples
get_percentile <- function(Exp.Data,
                           Meta.Data,
                           percentile){
  Exp.allGene <- get_AveExp(Gene = row.names(Exp.Data),
                            Exp.Data = Exp.Data,
                            Meta.Data = Meta.Data)
  N <- length(levels(Meta.Data$Type))

  if(length(percentile) == 0){
    percentile <- c(.5, .75, .90)
  }

  Np <- length(percentile)

  for (i in 1:N) {

    if(i==1){
      data <- data.frame(exp = Exp.allGene[, 1+i][which(Exp.allGene[, 1+i] > 0)])
      quantile <- quantile(data$exp, probs = percentile)
      quantile_1<-data.frame(t(matrix(quantile)))
      colnames(quantile_1)<-names(quantile)
    }else{
      data <- data.frame(exp = Exp.allGene[, 1+i][which(Exp.allGene[, 1+i] > 0)])
      quantile <- quantile(data$exp, probs = percentile)
      quantile_2<-data.frame(t(matrix(quantile)))
      colnames(quantile_2)<-names(quantile)

      quantile_1 <- rbind(quantile_1, quantile_2)
    }
  }
  row.names(quantile_1) <- levels(Meta.Data$Type)

  for (j in 1:Np) {
    if (j == 1){
      data1 <- data.frame(Type = row.names(quantile_1), Ave.Exp. = quantile_1[, j], Prob. = percentile[j])
    }else{
      data2 <- data.frame(Type = row.names(quantile_1), Ave.Exp. = quantile_1[, j], Prob. = percentile[j])
      data1 <- rbind(data1, data2)
    }
  }

  data1$Prob. <- as.factor(data1$Prob.)

  return(data1)
}

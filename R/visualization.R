#' plot heatmap of the gene expressions in clusters
#'
#' @param Gene exSigNet object (gene expressions and signaling strengths are optional)
#' @param Exp.Data expression matrix (rows: features, columns: cells)
#' @param Meta.Data meta data of the expression matrix, must have a variable "Type" containing all cell labels
#' @param cluster_rows boolean values determining if rows should be clustered or hclust object, the default value is "FALSE"
#' @param cluster_cols boolean values determining if columns should be clustered or hclust object, the default value is "FALSE"
#' @param cellwidth individual cell width in points, the default value is 15
#' @param cellheight individual cell height in points, the default value is 5
#' @param gaps_col similar to gaps_row, but for columns, the default value is "NULL"
#' @param gaps_row vector of row indices that show where to put gaps into heatmap. Used only if the rows are not clustered, the default value is "NULL"
#' @param fontsize_row fontsize for rownames, the default value is 10
#' @param fontsize_col fontsize for colnames, the default value is 10
#' @param angle_col angle of the column labels, right now one can choose only from few predefined options (0, 45, 90, 270 and 315), the default value is 0
#' @param legend_breaks vector of breakpoints for the legend
#' @param legend_labels vector of labels for the legend_breaks
#' @param annotation_row data frame that specifies the annotations shown on left side of the heatmap. Each row defines the features for a specific row. The rows in the data and in the annotation are matched using corresponding row names. Note that color schemes takes into account if variable is continuous or discrete
#' @param annotation_col similar to annotation_row, but for columns
#' @importFrom pheatmap pheatmap
#'
#' @return
#' @export
#'
#' @examples
get_heatmap <- function(Gene,
                        Exp.Data,
                        Meta.Data,
                        cluster_rows = FALSE,
                        cluster_cols = FALSE,
                        cellwidth = 15,
                        cellheight = 5,
                        gaps_col = NULL,
                        gaps_row = NULL,
                        fontsize_row = 10,
                        fontsize_col = 10,
                        angle_col = 0,
                        legend_breaks = NULL,
                        legend_labels = NULL,
                        annotation_row = NULL,
                        annotation_col = NULL){

  exS_exp <- get_AveExp(Gene = Gene,
                        Exp.Data = Exp.Data,
                        Meta.Data = Meta.Data)

  exS_exp.HP <- exS_exp
  exS_exp.HP <- exS_exp.HP[, -1]
  row.names(exS_exp.HP) <- exS_exp$Gene

  plot <- pheatmap(mat = exS_exp.HP,
                   cluster_rows = cluster_rows,
                   cluster_cols = cluster_cols,
                   cellwidth = cellwidth,
                   cellheight = cellheight,
                   gaps_col = NULL,
                   gaps_row = NULL,
                   fontsize_row = 10,
                   fontsize_col = 10,
                   angle_col = 90,
                   legend_breaks = NULL,
                   legend_labels = NULL,
                   annotation_row = NULL,
                   annotation_col = NULL)

  return(plot)
}

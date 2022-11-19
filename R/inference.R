#' infer ligand-target GRN from user-select target genes based on prior knowledge (exFINDER-DB)
#'
#' @param Target user-select target genes (a list of gene symbols)
#' @param DB the exFINDER-DB for this species (human, mouse, or zebrafish)
#'
#' @return
#' @export
#'
#' @examples
get_ltGRN <- function(Target,
                      DB){

  if (length(Target) == 0) {
    warning("Missing target genes!")
    return(NA_real_)}

  d.Target <- Target[duplicated(Target)]
  if (length(d.Target) != 0) {
    warning("Duplicated target genes found and ignored!")}

  Layer.TFT <- DB[[3]][DB[[3]]$to %in% Target, ]
  TF <- unique(DB[[3]][DB[[3]]$to %in% Target, ]$from)
  Layer.RTF <- DB[[2]][DB[[2]]$to %in% TF, ]
  R <- unique(DB[[2]][DB[[2]]$to %in% TF, ]$from)
  TF.2 <- unique(DB[[2]][DB[[2]]$to %in% TF, ]$to)
  Layer.LR <- DB[[1]][DB[[1]]$to %in% R, ]
  L <- unique(DB[[1]][DB[[1]]$to %in% R, ]$from)
  R.2 <- unique(DB[[1]][DB[[1]]$to %in% R, ]$to)

  if (length(L) == 0) {
    warning("No ligand inferred!")
    return(NA_real_)}

  L.2 <- unique(L)
  Layer1 <- DB[[1]][ (DB[[1]]$from %in% L.2) & (DB[[1]]$to %in% R.2), ]
  Layer1$from.Role <- "Ligand"
  Layer1$to.Role <- "Receptor"

  Layer2 <- DB[[2]][ (DB[[2]]$from %in% R.2) & (DB[[2]]$to %in% TF.2), ]
  Layer2$from.Role <- "Receptor"
  Layer2$to.Role <- "TF"

  Layer3 <- DB[[3]][(DB[[3]]$from %in% TF.2) & (DB[[3]]$to %in% Target), ]
  Layer3$from.Role <- "TF"
  Layer3$to.Role <- "Target"
  T.2 <- unique(Layer3$to)

  T.diff <- setdiff(Target, T.2)
  if (length(T.diff) == 0) {
    warning("The following target genes have been ignored: ", T.diff ,"\n")}

  Graph.Edge <- rbind(Layer1, Layer2, Layer3)
  rownames(Graph.Edge) <- 1:nrow(Graph.Edge)

  Node1 <- data.frame(Gene = L.2, Role = "Ligand")
  Node2 <- data.frame(Gene = R.2, Role = "Receptor")
  Node3 <- data.frame(Gene = TF.2, Role = "TF")
  Node4 <- data.frame(Gene = T.2, Role = "Target")
  Graph.Node <- rbind(Node1, Node2, Node3, Node4)

  NodeT <- c(L.2, R.2, TF.2, T.2)
  NodeT <- unique(NodeT)
  OverLap <- length(L.2) + length(R.2) + length(TF.2) + length(T.2) - length(NodeT)
  if (OverLap != 0) {
    warning("Overlap found!")}else {print("No overlap found!")}

  output = list()
  output[[1]] = Graph.Edge
  output[[2]] = Graph.Node
  names(output)<-c('Graph.Edge', 'Graph.Node')
  return(output)

}

#' infer the potential external signals based on ltGRN and expression data
#'
#' uses ltGRN, scRNA-seq data and meta data as inputs to (1) filter out the non-measured genes, (2) calculates average expressions of the ligands in every cluster, (3) infers potential external signals based on their expressions
#'
#' @param Graph ltGRN object (output of "get_ltGRN")
#' @param Exp.Data expression matrix (rows: features, columns: cells)
#' @param Meta.Data meta data of the expression matrix, must have a variable "Type" containing all cell labels
#' @param cutoff a threshold to filter the ligands. Negative cutoff value: selecting the ligands below this threshold; zero cutoff value: keep all ligands; otherwise selecting the ligands above this threshold
#' @param AG.L list of cell groups that forces exFINDER to filter the ligands based on these cell groups (optional)
#'
#' @return
#' @export
#'
#' @examples
get_potentialex <- function(Graph,
                            Exp.Data,
                            Meta.Data,
                            cutoff,
                            AG.L){

  Node <- Graph$Graph.Node
  Edge <- Graph$Graph.Edge

  Node.L <- Node$Gene[Node$Role == "Ligand"]
  Node.R <- Node$Gene[Node$Role == "Receptor"]
  Node.TF <- Node$Gene[Node$Role == "TF"]
  Node.T <- Node$Gene[Node$Role == "Target"]

  Node.TF <- intersect(Node.TF, row.names(Exp.Data))
  Node.R <- intersect(Node.R, row.names(Exp.Data))
  Node.L <- intersect(Node.L, row.names(Exp.Data))

  Edge.3 <- Edge[Edge$from.Role == "TF", ]
  Edge.3 <- Edge.3[ (Edge.3$from %in% Node.TF) & (Edge.3$to %in% Node.T), ]

  Edge.2 <- Edge[Edge$from.Role == "Receptor", ]
  Edge.2 <- Edge.2[ (Edge.2$from %in% Node.R) & (Edge.2$to %in% Node.TF), ]

  Edge.1 <- Edge[Edge$from.Role == "Ligand", ]
  Edge.1 <- Edge.1[ (Edge.1$from %in% Node.L) & (Edge.1$to %in% Node.R), ]

  L <- unique(Edge.1$from[Edge.1$from.Role == "Ligand"])

  Exp.L <- get_AveExp(Gene = L,
                      Exp.Data = Exp.Data,
                      Meta.Data = Meta.Data)
  if (length(AG.L)!=0){
    Exp.L <- Exp.L[, c('Gene', AG.L)]
    Exp.L2 <- get_RowMax(Exp.L)
  } else if (length(AG.L)==0){
    Exp.L2 <- get_RowMax(Exp.L)
  }

  if (cutoff < 0){
    L2 <- Exp.L2$Gene[Exp.L2$Max_AveExp < -cutoff]
  } else if (cutoff > 0){
    L2 <- Exp.L2$Gene[Exp.L2$Max_AveExp > cutoff]
  } else if (cutoff == 0){
    L2 <- Exp.L2$Gene
  }

  return(L2)
}

#' infer the external signal-target signaling network (exSigNet)
#'
#' uses scRNA-seq data and ltGRN to infer highly expressed TF and R, identify external signals, and infer the exSigNet
#'
#' @param Graph ltGRN object (output of "get_ltGRN")
#' @param Ligands potential external signals
#' @param Exp.Data expression matrix (rows: features, columns: cells)
#' @param Meta.Data meta data of the expression matrix, must have a variable "Type" containing all cell labels
#' @param par.cutoff thresholds to filter the ligands. Negative cutoff value: selecting the ligands below this threshold; zero cutoff value: keep all ligands; otherwise selecting the ligands above this threshold
#' @param AG.R list of cell groups that forces exFINDER to filter the receptors based on these cell groups (optional)
#' @param AG.TF list of cell groups that forces exFINDER to filter the TFs based on these cell groups (optional)
#'
#' @return
#' @export
#'
#' @examples
get_exSigNet <- function(Graph,
                         Ligands,
                         Exp.Data,
                         Meta.Data,
                         par.cutoff,
                         AG.R,
                         AG.TF){

  if (length(par.cutoff)==1){
    par.cutoff <- rep(par.cutoff, 2)
  } else if (length(par.cutoff)==0){
    warning("par.cutoff must have at least one cutoff value!")
    return(NA_real_)}

  Node <- Graph$Graph.Node
  Edge <- Graph$Graph.Edge

  Node.L <- Node$Gene[Node$Role == "Ligand"]
  Node.R <- Node$Gene[Node$Role == "Receptor"]
  Node.TF <- Node$Gene[Node$Role == "TF"]
  Node.T <- Node$Gene[Node$Role == "Target"]

  Node.TF <- intersect(Node.TF, row.names(Exp.Data))
  Node.R <- intersect(Node.R, row.names(Exp.Data))
  Node.L <- intersect(Node.L, row.names(Exp.Data))

  # filter TF
  Exp.TF <- get_AveExp(Gene = Node.TF,
                       Exp.Data = Exp.Data,
                       Meta.Data = Meta.Data)
  if (length(AG.TF)!=0){
    Exp.TF <- Exp.TF[, c('Gene', AG.TF)]
    Exp.TF <- get_RowMax(Exp.TF)
  } else if (length(AG.TF)==0){
    Exp.TF <- get_RowMax(Exp.TF)
  }

  if (par.cutoff[2] < 0){
    TF.2 <- Exp.TF$Gene[Exp.TF$Max_AveExp < -par.cutoff[2]]
  } else if (par.cutoff[2] > 0){
    TF.2 <- Exp.TF$Gene[Exp.TF$Max_AveExp > par.cutoff[2]]
  } else if (par.cutoff[2] == 0){
    TF.2 <- Exp.TF$Gene
  }

  if (length(TF.2)==0){
    warning("no such TF inferred!")
    return(NA_real_)}

  Edge.3 <- Edge[Edge$from.Role == "TF", ]
  Edge.3 <- Edge.3[ (Edge.3$from %in% TF.2) & (Edge.3$to %in% Node.T), ]

  # filter R
  Edge.2 <- Edge[Edge$from.Role == "Receptor", ]
  Edge.2 <- Edge.2[ (Edge.2$from %in% Node.R) & (Edge.2$to %in% TF.2), ]
  Node.R <- unique(Edge.2$from)

  Exp.R <- get_AveExp(Gene = Node.R,
                      Exp.Data = Exp.Data,
                      Meta.Data = Meta.Data)
  if (length(AG.R)!=0){
    Exp.R <- Exp.R[, c('Gene', AG.R)]
    Exp.R <- get_RowMax(Exp.R)
  } else if (length(AG.R)==0){
    Exp.R <- get_RowMax(Exp.R)
  }


  if (par.cutoff[1] < 0){
    R.2 <- Exp.R$Gene[Exp.R$Max_AveExp < -par.cutoff[1]]
  } else if (par.cutoff[1] > 0){
    R.2 <- Exp.R$Gene[Exp.R$Max_AveExp > par.cutoff[1]]
  } else if (par.cutoff[1] == 0){
    R.2 <- Exp.R$Gene
  }


  if (length(R.2)==0){
    warning("no such Receptor inferred!")
    return(NA_real_)}

  Edge.1 <- Edge[Edge$from.Role == "Ligand", ]
  Edge.1 <- Edge.1[ (Edge.1$from %in% Ligands) & (Edge.1$to %in% R.2), ]
  exS <- unique(Edge.1$from)

  Edge.1 <- Edge.1[ (Edge.1$from %in% exS) & (Edge.1$to %in% R.2), ]
  R.3 <- unique(Edge.1$to)
  Edge.2 <- Edge.2[ (Edge.2$from %in% R.3) & (Edge.2$to %in% TF.2), ]
  TF.3 <- unique(Edge.2$to)
  Edge.3 <- Edge.3[ (Edge.3$from %in% TF.3) & (Edge.3$to %in% Node.T), ]
  T.3 <- unique(Edge.3$to)

  Graph.Edge <- rbind(Edge.1, Edge.2, Edge.3)
  rownames(Graph.Edge) <- 1:nrow(Graph.Edge)

  Node1 <- data.frame(Gene = exS, Role = "Ligand")
  Node2 <- data.frame(Gene = R.3, Role = "Receptor")
  Node3 <- data.frame(Gene = TF.3, Role = "TF")
  Node4 <- data.frame(Gene = T.3, Role = "Target")
  Graph.Node <- rbind(Node1, Node2, Node3, Node4)

  NodeT <- c(exS, R.3, TF.3, T.3)
  NodeT <- unique(NodeT)
  OverLap <- length(exS) + length(R.3) + length(TF.3) + length(T.3) - length(NodeT)
  if (OverLap != 0) {
    warning("Overlap found!")}else {print("No overlap found!")}

  output = list()
  output[[1]] = Graph.Edge
  output[[2]] = Graph.Node
  names(output)<-c('Graph.Edge', 'Graph.Node')
  return(output)
}

#' fix the overlap of the external signal-target signaling network (exSigNet) if needed
#'
#' @param Graph exSigNet object (output of "get_exSigNet")
#'
#' @return
#' @export
#'
#' @examples
fix_OverLap <- function(Graph){

  Node <- Graph$Graph.Node
  Edge <- Graph$Graph.Edge

  Node.L <- unique(Node$Gene[Node$Role == "Ligand"])
  Node.R <- unique(Node$Gene[Node$Role == "Receptor"])
  Node.TF <- unique(Node$Gene[Node$Role == "TF"])
  Node.T <- unique(Node$Gene[Node$Role == "Target"])

  Edge.3 <- Edge[Edge$from.Role == "TF", ]
  Edge.3 <- Edge.3[ (Edge.3$from %in% Node.TF) & (Edge.3$to %in% Node.T), ]

  Edge.2 <- Edge[Edge$from.Role == "Receptor", ]
  Edge.2 <- Edge.2[ (Edge.2$from %in% Node.R) & (Edge.2$to %in% Node.TF), ]

  Edge.1 <- Edge[Edge$from.Role == "Ligand", ]
  Edge.1 <- Edge.1[ (Edge.1$from %in% Node.L) & (Edge.1$to %in% Node.R), ]

  common.L <- intersect(Node.L, c(Node.R, Node.TF, Node.T))
  if (length(common.L)!=0){
    Node.L <- setdiff(Node.L, common.L)
  }

  common.RT <- intersect(Node.R, Node.T)
  Node.R <- setdiff(Node.R, Node.T)
  common.TFT <- intersect(Node.TF, Node.T)
  Node.TF <- setdiff(Node.TF, Node.T)

  common.RTF <- intersect(Node.R, Node.TF)
  if (length(common.RTF)!=0){
    for (i in 1:length(common.RTF)) {
      node <- common.RTF[i]
      degree.R <- length(Edge.1$to[ (Edge.1$to == node) ]) + length(Edge.2$from[ (Edge.2$from == node) ])
      degree.TF <- length(Edge.2$to[ (Edge.2$to == node) ]) + length(Edge.3$from[ (Edge.3$from == node) ])
      if (degree.R < degree.TF){
        Node.R <- setdiff(Node.R, node)
      }else {Node.TF <- setdiff(Node.TF, node)}
    }
  }

  Edge.1 <- Edge.1[ (Edge.1$from %in% Node.L) & (Edge.1$to %in% Node.R) , ]
  Edge.2 <- Edge.2[ (Edge.2$from %in% Node.R) & (Edge.2$to %in% Node.TF), ]
  Edge.3 <- Edge.3[ (Edge.3$from %in% Node.TF) & (Edge.3$to %in% Node.T), ]

  Graph.Edge <- rbind(Edge.1, Edge.2, Edge.3)
  rownames(Graph.Edge) <- 1:nrow(Graph.Edge)

  Node1 <- data.frame(Gene = Node.L, Role = "Ligand")
  Node2 <- data.frame(Gene = Node.R, Role = "Receptor")
  Node3 <- data.frame(Gene = Node.TF, Role = "TF")
  Node4 <- data.frame(Gene = Node.T, Role = "Target")
  Graph.Node <- rbind(Node1, Node2, Node3, Node4)

  NodeT <- c(Node.L, Node.R, Node.TF, Node.T)
  NodeT <- unique(NodeT)
  OverLap <- length(Node.L) + length(Node.R) + length(Node.TF) + length(Node.T) - length(NodeT)
  if (OverLap != 0) {
    warning("Overlap found!")}else {print("No overlap found!")}

  output = list()
  output[[1]] = Graph.Edge
  output[[2]] = Graph.Node
  names(output)<-c('Graph.Edge', 'Graph.Node')
  return(output)
}

#' infer external signals (and exSigNet) that come from user-defined external cells
#'
#' @param Graph exSigNet object (output of "get_exSigNet")
#' @param Exp.ExData expression matrix of the external cells (rows: features, columns: cells)
#' @param Meta.ExData meta data of the expression matrix, must have a variable "Type" containing all cell labels
#' @param cutoff.ExData a threshold to filter the ligands (based on the Exp.ExData). Negative cutoff value: selecting the ligands below this threshold; zero cutoff value: keep all ligands; otherwise selecting the ligands above this threshold
#' @param AG.ExData list of cell groups (in the external cells) that forces exFINDER to filter the ligands based on these cell groups (optional)
#'
#' @return
#' @export
#'
#' @examples
get_exSigNetExData <- function(Graph,
                               Exp.ExData,
                               Meta.ExData,
                               cutoff.ExData,
                               AG.ExData){
  Node <- Graph$Graph.Node
  Edge <- Graph$Graph.Edge

  Node.L <- unique(Node$Gene[Node$Role == "Ligand"])
  Node.R <- unique(Node$Gene[Node$Role == "Receptor"])
  Node.TF <- unique(Node$Gene[Node$Role == "TF"])
  Node.T <- unique(Node$Gene[Node$Role == "Target"])

  Edge.3 <- Edge[Edge$from.Role == "TF", ]
  Edge.3 <- Edge.3[ (Edge.3$from %in% Node.TF) & (Edge.3$to %in% Node.T), ]

  Edge.2 <- Edge[Edge$from.Role == "Receptor", ]
  Edge.2 <- Edge.2[ (Edge.2$from %in% Node.R) & (Edge.2$to %in% Node.TF), ]

  Edge.1 <- Edge[Edge$from.Role == "Ligand", ]
  Edge.1 <- Edge.1[ (Edge.1$from %in% Node.L) & (Edge.1$to %in% Node.R), ]

  exS.Exp <- get_AveExp(Gene = Node.L,
                        Exp.Data = Exp.ExData,
                        Meta.Data = Meta.ExData)
  if (length(AG.ExData)!=0){
    exS.Exp <- exS.Exp[, c('Gene', AG.ExData)]
    exS.Exp <- get_RowMax(exS.Exp)
  } else if (length(AG.ExData)==0){
    exS.Exp <- get_RowMax(exS.Exp)
  }

  if (cutoff.ExData < 0){
    exS <- exS.Exp$Gene[exS.Exp$Max_AveExp < -cutoff.ExData]
  } else if (cutoff.ExData > 0){
    exS <- exS.Exp$Gene[exS.Exp$Max_AveExp > cutoff.ExData]
  } else if (cutoff.ExData == 0){
    exS <- exS.Exp$Gene
  }

  Edge.1 <- Edge.1[ (Edge.1$from %in% exS) & (Edge.1$to %in% Node.R) , ]
  R.3 <- unique(Edge.1$to)
  Edge.2 <- Edge.2[ (Edge.2$from %in% R.3) & (Edge.2$to %in% Node.TF), ]
  TF.3 <- unique(Edge.2$to)
  Edge.3 <- Edge.3[ (Edge.3$from %in% TF.3) & (Edge.3$to %in% Node.T), ]
  T.3 <- unique(Edge.3$to)

  Graph.Edge <- rbind(Edge.1, Edge.2, Edge.3)
  rownames(Graph.Edge) <- 1:nrow(Graph.Edge)

  Node1 <- data.frame(Gene = exS, Role = "Ligand")
  Node2 <- data.frame(Gene = R.3, Role = "Receptor")
  Node3 <- data.frame(Gene = TF.3, Role = "TF")
  Node4 <- data.frame(Gene = T.3, Role = "Target")
  Graph.Node <- rbind(Node1, Node2, Node3, Node4)

  NodeT <- c(exS, R.3, TF.3, T.3)
  NodeT <- unique(NodeT)
  OverLap <- length(exS) + length(R.3) + length(TF.3) + length(T.3) - length(NodeT)
  if (OverLap != 0) {
    warning("Overlap found!")}else {print("No overlap found!")}

  output = list()
  output[[1]] = Graph.Edge
  output[[2]] = Graph.Node
  names(output)<-c('Graph.Edge', 'Graph.Node')
  return(output)
}

#' find the receptors that only activated by external signals
#'
#' @param Graph exSigNet object (gene expressions and signaling strengths are optional)
#' @param Exp.Data expression matrix (rows: features, columns: cells)
#' @param DB the exFINDER-DB for this species (human, mouse, or zebrafish)
#'
#' @return
#' @export
#'
#' @examples
filterLR_exSigNet <- function(Graph,
                              Exp.Data,
                              DB){

  Node <- Graph$Graph.Node
  Edge <- Graph$Graph.Edge

  Node.L <- Node$Gene[Node$Role == "Ligand"]
  Node.R <- Node$Gene[Node$Role == "Receptor"]
  Node.TF <- Node$Gene[Node$Role == "TF"]
  Node.T <- Node$Gene[Node$Role == "Target"]

  n <- length(Node.R)
  R.filtered <- c()

  for (i in 1:n) {
    R.gene <- Node.R[i]
    L.gene <- unique(DB[[1]][DB[[1]]$to == R.gene, ]$from)
    L.gene <- intersect(L.gene, row.names(Exp.Data))
    if (length(setdiff(L.gene, Node.L)) == 0) {
      R.filtered <- c(R.filtered, R.gene)
    }
  }

  if (length(R.filtered)==0){
    warning("no such receptor found!")
    return(NA_real_)}

  Edge.1 <- Edge[Edge$from.Role == "Ligand", ]
  Edge.1 <- Edge.1[ (Edge.1$from %in% Node.L) & (Edge.1$to %in% R.filtered), ]
  exS <- unique(Edge.1$from)
  R.3 <- unique(R.filtered)

  Edge.2 <- Edge[Edge$from.Role == "Receptor", ]
  Edge.2 <- Edge.2[ (Edge.2$from %in% R.3) & (Edge.2$to %in% Node.TF), ]
  TF.3 <- unique(Edge.2$to)
  Edge.3 <- Edge[Edge$from.Role == "TF", ]
  Edge.3 <- Edge.3[ (Edge.3$from %in% TF.3) & (Edge.3$to %in% Node.T), ]
  T.3 <- unique(Edge.3$to)

  Graph.Edge <- rbind(Edge.1, Edge.2, Edge.3)
  rownames(Graph.Edge) <- 1:nrow(Graph.Edge)

  Node1 <- data.frame(Gene = exS, Role = "Ligand")
  Node2 <- data.frame(Gene = R.3, Role = "Receptor")
  Node3 <- data.frame(Gene = TF.3, Role = "TF")
  Node4 <- data.frame(Gene = T.3, Role = "Target")
  Graph.Node <- rbind(Node1, Node2, Node3, Node4)

  NodeT <- c(exS, R.3, TF.3, T.3)
  NodeT <- unique(NodeT)
  OverLap <- length(exS) + length(R.3) + length(TF.3) + length(T.3) - length(NodeT)
  if (OverLap != 0) {
    warning("Overlap found!")}else {print("No overlap found!")}

  output = list()
  output[[1]] = Graph.Edge
  output[[2]] = Graph.Node
  names(output)<-c('Graph.Edge', 'Graph.Node')
  return(output)
}

#' find the sub-exSigNet based on additional TF information
#'
#' @param Graph exSigNet object (gene expressions and signaling strengths are optional)
#' @param TF user-select TFs
#'
#' @return
#' @export
#'
#' @examples
get_subNet <- function(Graph,
                       TF){

  Node <- Graph$Graph.Node
  Edge <- Graph$Graph.Edge

  Node.L <- Node$Gene[Node$Role == "Ligand"]
  Node.R <- Node$Gene[Node$Role == "Receptor"]
  Node.TF <- Node$Gene[Node$Role == "TF"]
  Node.T <- Node$Gene[Node$Role == "Target"]

  if (length(TF) == 1){
    Node.TF2 <- Node$Gene[Node$Gene == TF]
  } else{
    Node.TF2 <- Node$Gene[Node$Gene %in% TF]
  }

  Edge.3 <- Edge[Edge$from.Role == "TF", ]
  Edge.3 <- Edge.3[ (Edge.3$from %in% Node.TF2) & (Edge.3$to %in% Node.T), ]
  Node.T3 <- unique(Edge.3$to)
  Node.TF3 <- unique(Edge.3$from)

  Edge.2 <- Edge[Edge$from.Role == "Receptor", ]
  Edge.2 <- Edge.2[ (Edge.2$from %in% Node.R) & (Edge.2$to %in% Node.TF3), ]
  Node.R3 <- unique(Edge.2$from)

  Edge.1 <- Edge[Edge$from.Role == "Ligand", ]
  Edge.1 <- Edge.1[ (Edge.1$from %in% Node.L) & (Edge.1$to %in% Node.R3), ]
  Node.L3 <- unique(Edge.1$from)

  if (length(Node.L3)==0 | length(Node.R3)==0 | length(Node.TF3)==0 | length(Node.T3)==0){
    warning("No such exSigNet inferred!")
    return(NA_real_)}

  E1 <- Edge[(Edge$from %in% Node.L3) & (Edge$to %in% Node.R3), ]
  E1$from.Role <- "Ligand"
  E1$to.Role <- "Receptor"
  E2 <- Edge[(Edge$from %in% Node.R3) & (Edge$to %in% Node.TF3), ]
  E2$from.Role <- "Receptor"
  E2$to.Role <- "TF"
  E3 <- Edge[(Edge$from %in% Node.TF3) & (Edge$to %in% Node.T3), ]
  E3$from.Role <- "TF"
  E3$to.Role <- "Target"
  Graph.Edge <- rbind(E1, E2, E3)
  rownames(Graph.Edge) <- 1:nrow(Graph.Edge)

  L1 <- Node[ Node$Gene %in% Node.L3, ]
  L2 <- Node[ Node$Gene %in% Node.R3, ]
  L3 <- Node[ Node$Gene %in% Node.TF3, ]
  L4 <- Node[ Node$Gene %in% Node.T3, ]
  Graph.Node <- rbind(L1, L2, L3, L4)

  output = list()
  output[[1]] = Graph.Edge
  output[[2]] = Graph.Node
  names(output)<-c('Graph.Edge', 'Graph.Node')

  return(output)
}

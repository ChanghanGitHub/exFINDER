#' show the details of a exSigNet (or ltGRN), number of nodes/edges, etc.
#'
#' @param Graph exSigNet object (output of "get_exSigNet")
#'
#' @return
#' @export
#'
#' @examples
get_NetworkDetail <- function(Graph){

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

  No.L <- length(Node.L)
  No.R <- length(Node.R)
  No.TF <- length(Node.TF)
  No.T <- length(Node.T)

  No.Edge1 <- length(Edge.1$from)
  No.Edge2 <- length(Edge.2$from)
  No.Edge3 <- length(Edge.3$from)

  cat("Number of ligands: ", No.L ,"\n")
  cat("Number of receptors: ", No.R ,"\n")
  cat("Number of tfs: ", No.TF ,"\n")
  cat("Number of targets: ", No.T ,"\n")

  cat("Number of L-R pairs: ", No.Edge1 ,"\n")
  cat("Number of R-TF pairs: ", No.Edge2 ,"\n")
  cat("Number of TF-T pairs: ", No.Edge3 ,"\n")
}

#' calculate the expression levels of Receptors, TFs, and Targets
#'
#' option: uses external dataset to calculate the expression levels of the external signals (otherwise they will be set to 1).
#'
#' @param Graph exSigNet object (output of "get_exSigNet")
#' @param Exp.Data expression matrix (rows: features, columns: cells)
#' @param Meta.Data meta data of the expression matrix, must have a variable "Type" containing all cell labels
#' @param AG.R list of cell groups that forces exFINDER to filter the receptors based on these cell groups (optional)
#' @param AG.TF list of cell groups that forces exFINDER to filter the TFs based on these cell groups (optional)
#' @param AG.T list of cell groups that forces exFINDER to filter the targets based on these cell groups (optional)
#' @param Exp.ExData expression matrix of the external cells (rows: features, columns: cells)
#' @param Meta.ExData meta data of the expression matrix, must have a variable "Type" containing all cell labels
#' @param AG.ExData list of cell groups (in the external cells) that forces exFINDER to filter the ligands based on these cell groups (optional)
#'
#' @return
#' @export
#'
#' @examples
get_NodeValue <- function(Graph,
                          Exp.Data,
                          Meta.Data,
                          AG.R,
                          AG.TF,
                          AG.T,
                          Exp.ExData,
                          Meta.ExData,
                          AG.ExData){

  Node <- Graph$Graph.Node

  Node.L <- Node$Gene[Node$Role == "Ligand"]
  Node.R <- Node$Gene[Node$Role == "Receptor"]
  Node.TF <- Node$Gene[Node$Role == "TF"]
  Node.T <- Node$Gene[Node$Role == "Target"]

  # calculate receptor
  R.Exp <- get_AveExp(Gene = Node.R,
                      Exp.Data = Exp.Data,
                      Meta.Data = Meta.Data)
  if (length(AG.R)!=0){
    R.Exp <- R.Exp[, c('Gene', AG.R)]
    R.Exp <- get_RowMax(R.Exp)
  } else if (length(AG.R)==0){
    R.Exp <- get_RowMax(R.Exp)
  }
  Value.R <- R.Exp$Max_AveExp

  # calculate TF
  TF.Exp <- get_AveExp(Gene = Node.TF,
                       Exp.Data = Exp.Data,
                       Meta.Data = Meta.Data)
  if (length(AG.TF)!=0){
    TF.Exp <- TF.Exp[, c('Gene', AG.TF)]
    TF.Exp <- get_RowMax(TF.Exp)
  } else if (length(AG.TF)==0){
    TF.Exp <- get_RowMax(TF.Exp)
  }
  Value.TF <- TF.Exp$Max_AveExp

  # calculate target
  T.Exp <- get_AveExp(Gene = Node.T,
                      Exp.Data = Exp.Data,
                      Meta.Data = Meta.Data)
  if (length(AG.T)!=0){
    T.Exp <- T.Exp[, c('Gene', AG.T)]
    T.Exp <- get_RowMax(T.Exp)
  } else if (length(AG.T)==0){
    T.Exp <- get_RowMax(T.Exp)
  }
  Value.T <- T.Exp$Max_AveExp

  # calculate ligand (external signal)
  if ( (is.null(Exp.ExData) == FALSE)&(is.null(Meta.ExData) == FALSE) ){

    L.Exp <- get_AveExp(Gene = Node.L,
                        Exp.Data = Exp.ExData,
                        Meta.Data = Meta.ExData)

    if (length(AG.ExData)!=0){
      L.Exp <- L.Exp[, c('Gene', AG.ExData)]
      L.Exp <- get_RowMax(L.Exp)
    } else if (length(AG.ExData)==0){
      L.Exp <- get_RowMax(L.Exp)
    }
    Value.L <- L.Exp$Max_AveExp
  }else { Value.L <- rep(1, length(Node.L)) }

  Value <- c(Value.L, Value.R, Value.TF, Value.T)

  Node$Size <- Value
  Graph$Graph.Node <- Node

  return(Graph)
}

#' predict the signaling strengths (edge weights)
#'
#' @param Graph exSigNet object with gene expressions (output of "get_NodeValue")
#' @param Kh parameter for the mass action law (default value is 2)
#'
#' @return
#' @export
#'
#' @examples
get_EdgeWeight <- function(Graph,
                           Kh){

  if (is.null(Kh)==TRUE){
    Kh <- 2
  }

  from.gene <- Graph$Graph.Edge$from
  loc.from <- match(from.gene, Graph$Graph.Node$Gene)
  from.size <- Graph$Graph.Node$Size[loc.from]

  to.gene <- Graph$Graph.Edge$to
  loc.to <- match(to.gene, Graph$Graph.Node$Gene)
  to.size <- Graph$Graph.Node$Size[loc.to]

  weight <- (from.size*to.size)/(Kh + (from.size*to.size))
  Graph$Graph.Edge$Weight <- weight

  return(Graph)
}

#' calculate the maximal flow between individual ligand and target
#'
#' @param Graph exSigNet object with gene expressions and signaling strengths
#' @importFrom igraph graph_from_data_frame graph.maxflow
#'
#' @return
#' @export
#'
#' @examples
get_LRMaxFlow <- function(Graph){

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

  N.L <- length(Node.L)
  N.T <- length(Node.T)
  Source <- c()
  Sink <- c()
  MaxFlow <- c()

  for (i in 1:N.L) {
    for (j in 1:N.T) {
      Source <- c( Source, Node.L[i])
      Sink <- c( Sink, Node.T[j])

      Edge.1 <- Edge[Edge$from.Role == "Ligand", ]
      Edge.1 <- Edge.1[ (Edge.1$from == Node.L[i]) & (Edge.1$to %in% Node.R), ]
      R.2 <- unique(Edge.1$to)

      Edge.3 <- Edge[Edge$from.Role == "TF", ]
      Edge.3 <- Edge.3[ (Edge.3$from %in% Node.TF) & (Edge.3$to %in% Node.T[j]), ]
      TF.2 <- unique(Edge.3$from)

      Edge.2 <- Edge[Edge$from.Role == "Receptor", ]
      Edge.2 <- Edge.2[ (Edge.2$from %in% R.2) & (Edge.2$to %in% TF.2), ]

      if ((length(Edge.1$from)==0)|(length(Edge.2$from)==0)|(length(Edge.3$from)==0) ){
        MaxFlow <- c(MaxFlow, 0)
      } else{
        R.3 <- unique(Edge.2$from)
        TF.3 <- unique(Edge.2$to)
        Edge.1 <- Edge.1[ (Edge.1$from == Node.L[i]) & (Edge.1$to %in% R.3), ]
        Edge.3 <- Edge.3[ (Edge.3$from %in% TF.3) & (Edge.3$to %in% Node.T[j]), ]
        Edge.flow <- rbind(Edge.1, Edge.2, Edge.3)
        edges <- data.frame(from = Edge.flow$from, to = Edge.flow$to, weight = Edge.flow$Weight)

        Node1 <- data.frame(Gene = Node.L[i], Role = "Ligand")
        Node2 <- data.frame(Gene = R.3, Role = "Receptor")
        Node3 <- data.frame(Gene = TF.3, Role = "TF")
        Node4 <- data.frame(Gene = Node.T[j], Role = "Target")
        nodes <- rbind(Node1, Node2, Node3, Node4)

        G <- graph_from_data_frame(d=edges, vertices=nodes, directed=TRUE)
        Max = graph.maxflow(G, Node.L[i], Node.T[j], E(G)$weight)$value
        MaxFlow <- c(MaxFlow, Max)
      }

    }
  }

  output <- data.frame(Source = Source, Sink = Sink, MaxFlow = MaxFlow)
  return(output)
}

#' calculate activation level (AcI) of a given exSigNet
#'
#' @param Graph exSigNet object with gene expressions and signaling strengths
#' @param Par parameter vector (1x4) for calculating the AcI, the default values are 2
#'
#' @return
#' @export
#'
#' @examples
get_AcI <- function(Graph,
                    Par){

  if (length(Par)==1){
    Par <- rep(Par, 4)
  } else if (length(Par)!=4){
    print("Par must be a scale or 1x4 vector!")
  }

  L.sum <- sum(Graph$Graph.Node$Size[Graph$Graph.Node$Role == "Ligand"])
  R.sum <- sum(Graph$Graph.Node$Size[Graph$Graph.Node$Role == "Receptor"])
  TF.sum <- sum(Graph$Graph.Node$Size[Graph$Graph.Node$Role == "TF"])
  T.sum <- sum(Graph$Graph.Node$Size[Graph$Graph.Node$Role == "Target"])
  G.size <- nrow(Graph$Graph.Node)*nrow(Graph$Graph.Edge)

  LR.flow <- sum(Graph$Graph.Edge$Weight[Graph$Graph.Edge$from.Role == "Ligand"])
  RTF.flow <- sum(Graph$Graph.Edge$Weight[Graph$Graph.Edge$from.Role == "Receptor"])
  TFT.flow <- sum(Graph$Graph.Edge$Weight[Graph$Graph.Edge$from.Role == "TF"])

  AcI <- (G.size/(Par[1]+G.size))* (LR.flow/(Par[2]+LR.flow))*(RTF.flow/(Par[3]+RTF.flow))*(TFT.flow/(Par[4]+TFT.flow))

  return(AcI)
}

#' find the functionally related pathways in the exSigNet
#'
#' @param Graph exSigNet object with gene expressions and signaling strengths
#' @param LR.DB user-select L-R database (here we use CellChatDB)
#' @param Par parameters for calculating AcIs, the default values are 2
#'
#' @return
#' @export
#'
#' @examples
check_Pathway <- function(Graph,
                          LR.DB,
                          Par){

  # LR.DB: must be CellChatDB
  LR.DB$path <- paste(LR.DB$ligand, LR.DB$receptor, sep="*")
  Node <- Graph$Graph.Node
  Edge <- Graph$Graph.Edge

  Edge.2 <- Edge[(Edge$from.Role == "Ligand") & (Edge$to.Role == "Receptor"), ]
  Path <- Edge.2[(Edge.2$path %in% LR.DB$path), ]
  Pathway.name <- LR.DB[(LR.DB$path %in% Path$path), ]$pathway_name
  Pathway.name <- unique(Pathway.name)
  L.Path <- unique(Path$from)
  R.Path <- unique(Path$to)

  if ( length(L.Path)==0 | length(R.Path)==0){
    warning("no pathway found!")
    return(NA_real_)}

  Pathway.AcI <- c()
  Pathway.Name <- list()
  Pathway.Graph <- list()

  for (i in 1:length(L.Path)) {
    L.PathNet <- L.Path[i]
    R.PathNet <- Path[Path$from == L.PathNet, ]$to
    TF.PathNet <- Edge[(Edge$from %in% R.PathNet) & (Edge$from.Role == "Receptor") & (Edge$to.Role == "TF"), ]$to
    T.PathNet <- Edge[(Edge$from %in% TF.PathNet) & (Edge$from.Role == "TF") & (Edge$to.Role == "Target"), ]$to

    Gene.PathNet <- c(L.PathNet, R.PathNet, TF.PathNet, T.PathNet)
    Node.PathNet <- Node[Node$Gene %in% Gene.PathNet, ]

    Edge.PathNet1 <- Edge[(Edge$from %in% L.PathNet) & (Edge$to %in% R.PathNet), ]
    Edge.PathNet2 <- Edge[(Edge$from %in% R.PathNet) & (Edge$to %in% TF.PathNet), ]
    Edge.PathNet3 <- Edge[(Edge$from %in% TF.PathNet) & (Edge$to %in% T.PathNet), ]
    Edge.PathNet <- rbind(Edge.PathNet1, Edge.PathNet2, Edge.PathNet3)

    Graph.PathNet = list()
    Graph.PathNet[[1]] = Edge.PathNet
    Graph.PathNet[[2]] = Node.PathNet
    names(Graph.PathNet)<-c('Graph.Edge', 'Graph.Node')

    Pathway.Graph[[i]] <- Graph.PathNet

    AcI <- get_AcI(Graph = Graph.PathNet,
                   Par = Par)
    Pathway.AcI <- c(Pathway.AcI, AcI)

    Pathway.loop2 <- c()
    for (j in 1:length(R.PathNet)) {
      Path.loop <- paste(L.PathNet, R.PathNet[j], sep="*")
      Pathway.loop <- unique(LR.DB[(LR.DB$path == Path.loop), ]$pathway_name)
      Pathway.loop2 <- c(Pathway.loop2, Pathway.loop)
    }
    Pathway.Name[[i]] <- Pathway.loop2

  }

  Pathway.Graph[[(length(L.Path)+1)]] <- Pathway.AcI
  Pathway.Graph[[(length(L.Path)+2)]] <- Pathway.Name

  return(Pathway.Graph)
}

#' prepare the matrix for k-mean clustering analysis
#'
#' @param Graph exSigNet object with gene expressions and signaling strengths
#'
#' @return
#' @export
#'
#' @examples
get_SimilarityMatrix <- function(Graph){

  Node <- Graph$Graph.Node
  Edge <- Graph$Graph.Edge
  Edge$path.2 <- paste(Edge$from, Edge$to, sep="_")

  Node.L <- Node$Gene[Node$Role == "Ligand"]
  Node.R <- Node$Gene[Node$Role == "Receptor"]
  Node.TF <- Node$Gene[Node$Role == "TF"]
  Node.T <- Node$Gene[Node$Role == "Target"]

  Graph.list <- list()
  names <- c()
  k <- 0

  for (i in 1:length(Node.L)) {
    for (j in 1:length(Node.T)) {
      Edge.11 <- Edge[Edge$from.Role == "Ligand", ]
      Edge.11 <- Edge.11[ (Edge.11$from == Node.L[i]) , ]
      R.1 <- unique(Edge.11$to)

      Edge.21 <- Edge[Edge$from.Role == "Receptor", ]
      Edge.21 <- Edge.21[ (Edge.21$from %in% R.1) , ]
      RTF.1 <- unique(Edge.21$path)

      Edge.32 <- Edge[Edge$from.Role == "TF", ]
      Edge.32 <- Edge.32[  (Edge.32$to == Node.T[j]), ]
      TF.2 <- unique(Edge.32$from)

      Edge.22 <- Edge[Edge$from.Role == "Receptor", ]
      Edge.22 <- Edge.22[ (Edge.22$to %in% TF.2) , ]
      RTF.2 <- unique(Edge.22$path)

      RTF.common <- intersect(RTF.1, RTF.2)

      if (length(RTF.common)==0){break}
      else{
        k <- 1 + k
        Edge.23 <- Edge[Edge$from.Role == "Receptor", ]
        Edge.23 <- Edge.23[ (Edge.23$path %in% RTF.common) , ]
        R.3 <- unique(Edge.23$from)
        TF.3 <- unique(Edge.23$to)
        L.3 <- Node.L[i]
        T.3 <- Node.T[j]

        E1 <- Edge[(Edge$from == Node.L[i]) & (Edge$to %in% R.3), ]
        E1$from.Role <- "Ligand"
        E1$to.Role <- "Receptor"
        E2 <- Edge[(Edge$from %in% R.3) & (Edge$to %in% TF.3), ]
        E2$from.Role <- "Receptor"
        E2$to.Role <- "TF"
        E3 <- Edge[(Edge$from %in% TF.3) & (Edge$to == Node.T[j]), ]
        E3$from.Role <- "TF"
        E3$to.Role <- "Target"
        Graph.Edge <- rbind(E1, E2, E3)
        rownames(Graph.Edge) <- 1:nrow(Graph.Edge)

        L1 <- Node[ Node$Gene %in% L.3, ]
        L2 <- Node[ Node$Gene %in% R.3, ]
        L3 <- Node[ Node$Gene %in% TF.3, ]
        L4 <- Node[ Node$Gene %in% T.3, ]
        Graph.Node <- rbind(L1, L2, L3, L4)

        output = list()
        output[[1]] = Graph.Edge
        output[[2]] = Graph.Node
        names(output)<-c('Graph.Edge', 'Graph.Node')
        Graph.list[[k]] <- output

        names <-c(paste0(Node.L[i],"_",Node.T[j]), names)

      }
    }
  }

  names(Graph.list) <- names

  df.1 <- data.frame(Full_network = Node$Size)
  rownames(df.1) <- Node$Gene
  df.2 <- data.frame(Full_network = Edge$Weight)
  rownames(df.2) <- Edge$path.2
  df <- rbind(df.1, df.2)
  df <- data.frame(t(df))

  for (i in 1:length(Graph.list)) {
    Graph.list[[i]]$Graph.Edge$path.2 <- paste(Graph.list[[i]]$Graph.Edge$from, Graph.list[[i]]$Graph.Edge$to, sep="_")
    loc.1 <- match(Graph.list[[i]]$Graph.Node$Gene, colnames(df))
    loc.2 <- match(Graph.list[[i]]$Graph.Edge$path.2, colnames(df))

    df.add <- integer(length(df))
    df.add[loc.1] <- Graph.list[[i]]$Graph.Node$Size
    df.add[loc.2] <- Graph.list[[i]]$Graph.Edge$Weight
    df.add <- data.frame(t(df.add))
    colnames(df.add) <- colnames(df)
    rownames(df.add) <- names(Graph.list)[i]
    df <- rbind(df, df.add)
  }

  df <- df[-1, ]

  return(df)
}

#' perform enrichment analysis of a exSigNet
#'
#' @param Graph exSigNet object with gene expressions and signaling strengths
#' @param OrgDb annotation db
#' @param Number.of.Terms output GO terms based on the p-values
#' @param keyType keytype of input gene, default is "ENTREZID"
#' @param ont One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
#' @param pvalueCutoff adjusted pvalue cutoff on enrichment tests to report, default is 0.05
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param universe background genes. If missing, the all genes listed in the database (eg TERM2GENE table) will be used as background.
#' @param qvalueCutoff qvalue cutoff on enrichment tests to report as significant. Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues and iii) qvalueCutoff on qvalues to be reported. default is 0.2
#' @param minGSSize minimal size of genes annotated by Ontology term for testing, default is 10
#' @param maxGSSize maximal size of genes annotated for testing, default is 500
#' @param readable whether mapping gene ID to gene Name, default is "FALSE"
#' @param pool If ont='ALL', whether pool 3 GO sub-ontologies, default is "FALSE"
#' @importFrom clusterProfiler bitr enrichGO
#' @import org.Mm.eg.db
#' @import org.Dr.eg.db
#' @import org.Hs.eg.db
#'
#' @return
#' @export
#'
#' @examples
get_EnrichAnalysis <- function(Graph,
                               OrgDb,
                               Number.of.Terms,
                               keyType = "ENTREZID",
                               ont = "MF",
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH",
                               universe,
                               qvalueCutoff = 0.2,
                               minGSSize = 10,
                               maxGSSize = 500,
                               readable = FALSE,
                               pool = FALSE){

  n <- Number.of.Terms
  Gene <- Graph$Graph.Node$Gene
  Gene.entrezid <- bitr(Gene, fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = OrgDb,
                        drop = FALSE)

  Gene.entrezid = Gene.entrezid[!duplicated(Gene.entrezid$ENTREZID),]
  Gene.entrezid <- Gene.entrezid[complete.cases(Gene.entrezid$ENTREZID),]

  go_enrich <- enrichGO(gene = Gene.entrezid$ENTREZID,
                        OrgDb = OrgDb,
                        keyType = keyType,
                        ont = ont,
                        pvalueCutoff = pvalueCutoff,
                        pAdjustMethod = pAdjustMethod,
                        universe,
                        qvalueCutoff = qvalueCutoff,
                        minGSSize = minGSSize,
                        maxGSSize = maxGSSize,
                        readable = readable,
                        pool = pool)

  if (length(go_enrich@result$ID)<n){
    warning("Not enough GO terms are inferred!")
    return(NA_real_)}

  GO.ID <- go_enrich@result$ID[1:n]
  GO.Name <- go_enrich@result$Description[1:n]

  Node.Pct <- c()
  Exp.Pct <- c()

  for (i in 1:length(GO.ID)) {

    GO.Gene <- data.frame(go_enrich@geneSets[GO.ID[i]])
    colnames(GO.Gene) <- "GO.Gene"
    loc <- match(GO.Gene$GO.Gene, Gene.entrezid$ENTREZID)
    GO.GeneSym <- Gene.entrezid$SYMBOL[loc]
    GO.GeneSym <- GO.GeneSym[complete.cases(GO.GeneSym)]

    Index.1 <- length(GO.GeneSym)/length(Graph$Graph.Node$Gene)
    Index.2 <- sum(Graph$Graph.Node$Size[Graph$Graph.Node$Gene %in% GO.GeneSym])/sum(Graph$Graph.Node$Size)

    Exp.Pct <- c(Index.2, Exp.Pct)
    Node.Pct <- c(Index.1, Node.Pct)
  }

  df.GOIndex <- data.frame(Node.Pct = Node.Pct, Exp.Pct = Exp.Pct)
  rownames(df.GOIndex) <- GO.ID
  df.GOIndex$GO.Name <- GO.Name

  output <- list()
  output[[1]] = go_enrich
  output[[2]] = df.GOIndex
  names(output)<-c('GO analysis results', 'Evaluation')

  return(output)
}

#' calculate the total outflows(inflows) of ligands(targets)
#'
#' @param df dataframe, output of "get_LRMaxFlow"
#' @importFrom stats aggregate
#'
#' @return
#' @export
#'
#' @examples
compare_MaxFlow <- function(df){

  df.exSignals <- aggregate(df$MaxFlow, by=list(Gene = df$Source), sum)
  df.exSignals <- data.frame(Gene = df.exSignals$Gene, TotalFlow = df.exSignals$x, Role = "Ligand")

  df.targets <- aggregate(df$MaxFlow, by=list(Gene = df$Sink), sum)
  df.targets <- data.frame(Gene = df.targets$Gene, TotalFlow = df.targets$x, Role = "Target")

  df.TotalFlow <- rbind(df.exSignals, df.targets)
  return(df.TotalFlow)
}

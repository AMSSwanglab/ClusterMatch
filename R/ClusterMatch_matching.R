#' Matching clusters of two datasets
#' @param data the correlation between clusters of two datasets
#' @param top the maximum number of matches for each cluster
#' @param threshold the minimum threshold for matching
#' @return matching matrix

Match<-function(data,top=3,threshold=1.96){
  match_matrix<-data
  for (i in 1:nrow(data)) {
    for (j in 1:ncol(data)) {
      top_row<-data[i,order(data[i,],decreasing = TRUE)[1:top]]
      top_col<-data[order(data[,j],decreasing = TRUE)[1:top],j]
      if (data[i,j] %in% top_row & data[i,j] %in% top_col & data[i,j]>threshold) {
        match_matrix[i,j]<-1
      }
      else
        match_matrix[i,j]<-0
    }

  }

  for (i in 1:nrow(data)) {
    for (j in 1:ncol(data)) {
      top_row<-data[i,order(data[i,],decreasing = TRUE)[1]]
      top_col<-data[order(data[,j],decreasing = TRUE)[1],j]
      if (data[i,j] %in% top_row & data[i,j] %in% top_col & data[i,j]>threshold) {
        match_matrix[i,j]<-2
      }
    }

  }

  return(match_matrix)

}


#' Calculating the matching matrix between two datasets at the optimal resolutions
#' @param reference_matrix the gene expression matrix of dataset 1
#' @param query_matrix the gene expression matrix of dataset 2
#' @param ref.res the optimal clustering resolution of dataset 1
#' @param que.res the optimal clustering resolution of dataset 2
#' @param ref.norm whether dataset 1 is normalized
#' @param que.norm whether dataset 2 is normalized
#' @param dim.pca number of PCs to calculate (50 by default)
#' @param dim.cca number of canonical vectors to calculate (30 by default)
#' @param min.cluster.cells the minimum number of cells in a cluster
#' @param random_PCC the minimum value of random noise
#' @param min.pct Seurat FindAllMarkers min.pct
#' @param logfc.threshold Seurat FindAllMarkers logfc.threshold
#' @return a list containing the following elements:
#' - D1_reference: SeuratObject of dataset 1
#' - D2_query: SeuratObject of dataset 2
#' - Matching_matrix: matching matrix, 2 is mutually most correlated clusters, 1 is mutually correlated clusters
#' - strongest_matching_marker_top10: top 10 common marker genes of mutually most correlated clusters
#' - com_marker_list: common marker genes of mutually most correlated clusters
#' - Epsilon: random noise
#' - CCA_cor: PCC between clusters based on CCA representation
#' - Marker_cor：Jaccard coefficient between clusters based on CCA representation
#' - ClusterMatch_cor = ClusterMatch correlation between clusters
#' @importFrom Seurat CreateSeuratObject NormalizeData ScaleData RunPCA FindNeighbors FindClusters RunCCA L2Dim SplitObject FindAllMarkers Idents
#' @importFrom dplyr select left_join
#' @importFrom stringr str_c
#' @export

ClusterMatch_matching <- function(reference_matrix, query_matrix, ref.res, que.res, ref.norm = TRUE, que.norm = TRUE,
                                  dim.pca = 50, dim.cca = 30, min.cluster.cells = 20, random_PCC = 0, min.pct = 0.25, logfc.threshold = 0.25){
  #object and clustering
  reference <- Seurat::CreateSeuratObject(counts = reference_matrix)
  if (ref.norm == TRUE) {
    reference <- Seurat::NormalizeData(reference)
  }
  reference <- Seurat::FindVariableFeatures(reference, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  reference <- Seurat::ScaleData(reference, verbose = FALSE)
  reference <- Seurat::RunPCA(reference, features = Seurat::VariableFeatures(object = reference), verbose = FALSE)
  reference <- Seurat::FindNeighbors(reference, dims = 1:dim.pca,verbose = FALSE)
  reference <- Seurat::FindClusters(reference, resolution = ref.res)
  reference@meta.data$name <- "D1_"
  reference$seurat_clusters <- stringr::str_c(reference$name, reference$seurat_clusters)
  reference@meta.data <- dplyr::select(reference@meta.data, -c("name"))

  query <- Seurat::CreateSeuratObject(counts = query_matrix)
  if (que.norm == TRUE) {
    query <- Seurat::NormalizeData(query)
  }
  query <- Seurat::FindVariableFeatures(query, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  query <- Seurat::ScaleData(query, verbose = FALSE)
  query <- Seurat::RunPCA(query, features = Seurat::VariableFeatures(object = query), verbose = FALSE)
  query <- Seurat::FindNeighbors(query, dims = 1:dim.pca, verbose = FALSE)
  query <- Seurat::FindClusters(query, resolution = que.res)
  query@meta.data$name <- "D2_"
  query$seurat_clusters <- stringr::str_c(query$name,query$seurat_clusters)
  query@meta.data <- dplyr::select(query@meta.data, -c("name"))

  #CCA embedding
  CCA <- Seurat::RunCCA(object1 = reference, object2 = query, num.cc = dim.cca)
  L2CCA <- Seurat::L2Dim(CCA, reduction = "cca")
  embedding <- L2CCA@reductions[["cca.l2"]]@cell.embeddings
  embedding <- data.frame(cbind(cells=rownames(embedding), embedding))

  #CCA cluster representation
  #reference
  reference_cluster <- Seurat::SplitObject(reference, split.by = "seurat_clusters")
  reference_means <- data.frame()
  reference_cell_cluster_embedding <- list()
  for (i in 1:length(reference_cluster)) {
    cluster_name <- names(reference_cluster)[i]
    cluster_cell <- data.frame(cells = rownames(reference_cluster[[i]]@meta.data))
    cluster_embedding <- dplyr::left_join(cluster_cell, embedding, by="cells")
    rownames(cluster_embedding) <- cluster_embedding[, 1]
    cluster_embedding <- cluster_embedding[, -1]
    reference_cell_cluster_embedding[[i]] <- as.data.frame(lapply(cluster_embedding, as.numeric))
    rownames(reference_cell_cluster_embedding[[i]]) <- cluster_cell$cells
    reference_cell_cluster_embedding[[i]] <- t(reference_cell_cluster_embedding[[i]])
    names(reference_cell_cluster_embedding)[i] <- cluster_name
    cluster_embedding <- as.data.frame(lapply(cluster_embedding, as.numeric))
    cluster_vec <- colMeans(cluster_embedding)
    cluster_vec <- data.frame(cluster_vec)
    colnames(cluster_vec)[1] <- cluster_name
    if (i==1) {
      reference_means <- cluster_vec
    }
    else
      reference_means <- cbind(reference_means, cluster_vec)
  }
  #query
  query_cluster <- Seurat::SplitObject(query, split.by = "seurat_clusters")
  query_means <- data.frame()
  query_cell_cluster_embedding <- list()
  for (i in 1:length(query_cluster)) {
    cluster_name <- names(query_cluster)[i]
    cluster_cell <- data.frame(cells = rownames(query_cluster[[i]]@meta.data))
    cluster_embedding <- dplyr::left_join(cluster_cell, embedding, by="cells")
    rownames(cluster_embedding) <- cluster_embedding[, 1]
    cluster_embedding <- cluster_embedding[, -1]
    query_cell_cluster_embedding[[i]] <- as.data.frame(lapply(cluster_embedding, as.numeric))
    rownames(query_cell_cluster_embedding[[i]]) <- cluster_cell$cells
    query_cell_cluster_embedding[[i]] <- t(query_cell_cluster_embedding[[i]])
    names(query_cell_cluster_embedding)[i] <- cluster_name
    cluster_embedding <- as.data.frame(lapply(cluster_embedding, as.numeric))
    cluster_vec <- colMeans(cluster_embedding)
    cluster_vec <- data.frame(cluster_vec)
    colnames(cluster_vec)[1] <- cluster_name
    if (i==1) {
      query_means<-cluster_vec

    }
    else
      query_means<-cbind(query_means, cluster_vec)

  }

  #marker cluster representation
  #reference
  Seurat::Idents(object = reference) <- reference@meta.data$seurat_clusters
  reference_marker <- Seurat::FindAllMarkers(reference, only.pos = TRUE, min.pct = min.pct, logfc.threshold = logfc.threshold)
  reference_marker_cluster <- split(reference_marker, reference_marker$cluster)
  #query
  Seurat::Idents(object = query) <- query@meta.data$seurat_clusters
  query_marker <- Seurat::FindAllMarkers(query, only.pos = TRUE, min.pct = min.pct, logfc.threshold = logfc.threshold)
  query_marker_cluster <- split(query_marker,query_marker$cluster)

  #correlation of clusters
  #CCA
  CCA_PCC <- cor(query_means,reference_means)
  CCA_PCC[CCA_PCC<0] <- 0
  #marker
  jaccard <- function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
  }

  marker_jaccard<-data.frame()
  for (i in 1:length(query_marker_cluster)) {
    for (j in 1:length(reference_marker_cluster)) {
      marker_jaccard[i,j] <- jaccard(query_marker_cluster[[i]]$gene,reference_marker_cluster[[j]]$gene)

    }

  }
  rownames(marker_jaccard)<-names(query_marker_cluster)
  colnames(marker_jaccard)<-names(reference_marker_cluster)
  marker_jaccard[marker_jaccard=='NaN']<-0

  #Epsilon and beta
  reference_all_means <- embedding[1:ncol(reference_matrix), ]
  reference_all_means <- reference_all_means[, -1]
  reference_all_means <- as.data.frame(lapply(reference_all_means, as.numeric))
  reference_all_means_vec <- colMeans(reference_all_means)
  query_all_means <- embedding[(ncol(reference_matrix)+1):nrow(embedding), ]
  query_all_means <- query_all_means[, -1]
  query_all_means <- as.data.frame(lapply(query_all_means, as.numeric))
  query_all_means_vec <- colMeans(query_all_means)
  Epsilon <- cor(reference_all_means_vec, query_all_means_vec)
  beta <- max(CCA_PCC)/max(marker_jaccard)
  Epsilon <- max(Epsilon, random_PCC)
  CCA_marker_similarity <- CCA_PCC + beta*marker_jaccard
  CCA_marker_similarity[CCA_marker_similarity < Epsilon] <- 0
  match_matrix <- Match(data = CCA_marker_similarity, top = 3, threshold = Epsilon)
  CCA_PCC_match_matrix<-Match(data = CCA_PCC,top = 1,threshold = Epsilon)

  #strongest matching cluster marker
  overlap<-data.frame()
  flag=1
  for (i in 1:nrow(match_matrix)) {
    for (j in 1:ncol(match_matrix)) {
      if (match_matrix[i,j]==2) {
        overlap[flag,1] <- colnames(match_matrix)[j]
        overlap[flag,2] <- rownames(match_matrix)[i]
        flag<-flag+1


      }


    }

  }
  overlap$name <- c('/')
  overlap$name <- stringr::str_c(overlap$V1, overlap$name)
  overlap$name <- stringr::str_c(overlap$name, overlap$V2)
  reference_marker_strongest_matching_list <- subset(reference_marker_cluster, names(reference_marker_cluster) %in% overlap$V1)
  reference_marker_strongest_matching_list <- reference_marker_strongest_matching_list[overlap$V1]
  query_marker_strongest_matching_list <- subset(query_marker_cluster,names(query_marker_cluster) %in% overlap$V2)
  query_marker_strongest_matching_list <- query_marker_strongest_matching_list[overlap$V2]
  com_marker <- data.frame()
  com_marker_list <- list()
  com_marker_top_10 <- data.frame()
  for (i in 1:nrow(overlap)) {
    com_marker_list[[i]] <- merge(reference_marker_strongest_matching_list[[i]], query_marker_strongest_matching_list[[i]], by = "gene" )
    com_marker_list[[i]]$avg_log2FC <- (com_marker_list[[i]]$avg_log2FC.x+com_marker_list[[i]]$avg_log2FC.y)/2
    com_marker_list[[i]] <- com_marker_list[[i]][order(com_marker_list[[i]]$avg_log2FC, decreasing = TRUE), ]
    if (i==1) {
      com_marker_top_10 <- data.frame(com_marker_list[[i]]$gene[1:10])
    }
    else
      com_marker_top_10 <- cbind(com_marker_top_10, data.frame(com_marker_list[[i]]$gene[1:10]))
  }
  names(com_marker_list) <- overlap$name
  colnames(com_marker_top_10) <- overlap$name

  #output
  ClusterMatch_matching <- list(D1_reference = reference, D2_query = query, Matching_matrix = match_matrix,
                                strongest_matching_marker_top10 = com_marker_top_10, com_marker_list = com_marker_list, Epsilon = Epsilon, CCA_cor = CCA_PCC, Marker_cor = marker_jaccard, ClusterMatch_cor = CCA_marker_similarity )
}

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


#' Integrate scRNA-seq datasets
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
#' @param distance_diff regulate the different clusters distance
#' @param distance_same regulate the same clusters distance
#' @return a list containing the following elements:
#' - D1_reference: SeuratObject of dataset 1
#' - D2_query: SeuratObject of dataset 2
#' - Matching_matrix: matching matrix, 2 is mutually most correlated clusters, 1 is mutually correlated clusters
#' - cell_embedding: cell embeddings
#' @importFrom Seurat CreateSeuratObject NormalizeData ScaleData RunPCA FindNeighbors FindClusters RunCCA L2Dim SplitObject FindAllMarkers Idents
#' @importFrom dplyr select left_join
#' @importFrom stringr str_c
#' @export

ClusterMatch_integration <- function(reference_matrix, query_matrix, ref.res, que.res, ref.norm = TRUE, que.norm = TRUE,
                                     dim.pca = 50, dim.cca = 30, min.cluster.cells = 20, random_PCC = 0, distance_diff = 0, distance_same = 0){
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
  reference_marker <- Seurat::FindAllMarkers(reference, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  reference_marker_cluster <- split(reference_marker, reference_marker$cluster)
  #query
  Seurat::Idents(object = query) <- query@meta.data$seurat_clusters
  query_marker <- Seurat::FindAllMarkers(query, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
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
  rownames(marker_jaccard) <- names(query_marker_cluster)
  colnames(marker_jaccard) <- names(reference_marker_cluster)
  marker_jaccard[marker_jaccard=='NaN'] <- 0

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
  match_matrix1 <- match_matrix
  match_matrix[match_matrix == 2] <- 1
  match_matrix <- t(as.matrix(match_matrix))

  reference_adj_matrix <- matrix(0, nrow = nrow(match_matrix), ncol = nrow(reference@meta.data), dimnames = list(c(rownames(match_matrix)), c(rownames(reference@meta.data))))
  for (i in 1:ncol(reference_adj_matrix)) {
    cluster_id <- reference@meta.data$seurat_clusters[i]
    j <- which(rownames(match_matrix) == cluster_id)
    reference_adj_matrix[j, i] <- 1

  }
  query_adj_matrix <- matrix(0, nrow = ncol(match_matrix), ncol = nrow(query@meta.data), dimnames = list(c(colnames(match_matrix)), c(rownames(query@meta.data))))
  for (i in 1:ncol(query_adj_matrix)) {
    cluster_id <- query@meta.data$seurat_clusters[i]
    j <- which(colnames(match_matrix) == cluster_id)
    query_adj_matrix[j, i] <- 1

  }
  reference_adj_matrix <- t(reference_adj_matrix)
  query_adj_matrix <- t(query_adj_matrix)
  reference_cell_embedding <- embedding[1:nrow(reference@meta.data),]
  reference_cell_embedding <- reference_cell_embedding[,-1]
  reference_cell_embedding <- apply(as.matrix(reference_cell_embedding), 2, as.numeric)
  rownames(reference_cell_embedding) <- rownames(embedding)[1:nrow(reference@meta.data)]
  query_cell_embedding <- embedding[(nrow(reference@meta.data)+1):nrow(embedding),]
  query_cell_embedding <- query_cell_embedding[,-1]
  query_cell_embedding <- apply(as.matrix(query_cell_embedding), 2, as.numeric)
  rownames(query_cell_embedding) <- rownames(embedding)[(nrow(reference@meta.data)+1):nrow(embedding)]
  ones_matrix <- matrix(1, nrow = nrow(match_matrix), ncol = ncol(match_matrix) )
  reference_means <- as.matrix(reference_means)
  query_means <- as.matrix(query_means)

  #cluster correction
  reference_means1 <- reference_means + query_means %*% t(distance_diff * (match_matrix-ones_matrix) )
  query_means1 <- query_means + reference_means %*% (distance_diff * (match_matrix-ones_matrix) )
  reference_means2 <- reference_means1 + query_means1 %*% t(distance_same * match_matrix)
  query_means2 <- query_means1 + reference_means1 %*% (distance_same * match_matrix)

  reference_cell_embedding1 <- reference_cell_embedding + reference_adj_matrix %*% t(reference_means2 - reference_means)
  query_cell_embedding1 <- query_cell_embedding + query_adj_matrix %*% t(query_means2 - query_means)
  L2norm<-function(data){
    norm<-c()
    L2Data<-data
    z<-t(data) %*% data
    for (i in 1:ncol(data)) {
      norm[i]<-sqrt(z[i,i])

    }
    for (i in 1:ncol(data)) {
      L2Data[,i]<-data[,i]/norm[i]

    }
    rownames(L2Data)<-rownames(data)
    colnames(L2Data)<-colnames(data)
    return(L2Data)

  }
  cluster_means <- cbind(reference_means2, query_means2)
  L2_cluster_means <- L2norm(cluster_means)
  cell_embedding <- rbind(reference_cell_embedding1,query_cell_embedding1)
  L2_cell_embedding <- t(L2norm(t(cell_embedding)))

  #cell correct embedding
  cell_correct_embedding <- L2_cell_embedding %*% L2_cluster_means

  #output
  ClusterMatch_integration <- list(D1_reference = reference, D2_query = query, Matching_matrix = match_matrix1, cell_embedding = cell_correct_embedding)
  return(ClusterMatch_integration)
}

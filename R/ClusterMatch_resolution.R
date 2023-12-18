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


#' Calculating the optimal clustering resolutions for two datasets.
#' @param reference_matrix the gene expression matrix of dataset 1
#' @param query_matrix the gene expression matrix of dataset 2
#' @param ref.norm whether dataset 1 is normalized
#' @param que.norm whether dataset 2 is normalized
#' @param start.resolution starting resolution
#' @param end.resolution ending resolution
#' @param step step size between starting resolution and ending resolution
#' @param dim.pca number of PCs to calculate (50 by default)
#' @param dim.cca number of canonical vectors to calculate (30 by default)
#' @param min.cluster.cells the minimum number of cells in a cluster
#' @return a list containing the following elements:
#' - Resolution_cor: The average PCC of mutually most correlated clusters
#' - D1_ref_res: Optimal resolution for dataset 1
#' - D2_que_res: Optimal resolution for dataset 2
#' @importFrom Seurat CreateSeuratObject NormalizeData ScaleData RunPCA FindNeighbors FindClusters RunCCA L2Dim
#' @importFrom dplyr select
#' @importFrom stringr str_c
#' @export

ClusterMatch_resolution <- function(reference_matrix, query_matrix, ref.norm = TRUE, que.norm = TRUE, start.resolution = 0.1, end.resolution = 5, step = 0.1,
                                    dim.pca = 50, dim.cca = 30, min.cluster.cells = 20){
  #object
  reference <- Seurat::CreateSeuratObject(counts = reference_matrix)
  if (ref.norm == TRUE) {
    reference <- Seurat::NormalizeData(reference)
  }
  reference <- Seurat::FindVariableFeatures(reference, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  reference <- Seurat::ScaleData(reference,verbose = FALSE)
  reference <- Seurat::RunPCA(reference, features = Seurat::VariableFeatures(object = reference), verbose = FALSE)
  reference <- Seurat::FindNeighbors(reference, dims = 1:dim.pca, verbose = FALSE)

  query<-Seurat::CreateSeuratObject(counts = query_matrix)
  if (que.norm == TRUE) {
    query <- Seurat::NormalizeData(query)
  }
  query <- Seurat::FindVariableFeatures(query, selection.method = "vst", nfeatures = 2000,
                                        verbose = FALSE)
  query <- Seurat::ScaleData(query,verbose = FALSE)
  query <- Seurat::RunPCA(query, features = Seurat::VariableFeatures(object = query),verbose = FALSE)
  query <- Seurat::FindNeighbors(query, dims = 1:dim.pca,verbose = FALSE)

  #Louvain clustering
  for (i in seq(from=start.resolution, to=end.resolution, by = step)) {
    reference <- Seurat::FindClusters(reference, resolution = i)
    reference_cluster_cells <- split(reference@meta.data, reference@meta.data$seurat_clusters)
    num.cells <- c()
    min.cells <- c()
    for (j in 1:length(reference_cluster_cells)) {
      num.cells[j] <- nrow(reference_cluster_cells[[j]])
    }
    min.cells <- min(num.cells)
    if(min.cells < min.cluster.cells)
    {
      reference@meta.data <- reference@meta.data[,-ncol(reference@meta.data)]
      reference@meta.data$seurat_clusters <- reference@meta.data[,ncol(reference@meta.data)]
      break
    }
  }
  reference@meta.data <- dplyr::select(reference@meta.data, -c("seurat_clusters"))
  reference@meta.data$name <- "D1_"
  for (i in 4:(ncol(reference@meta.data)-1)) {
    reference@meta.data[, i] <- stringr::str_c(reference$name, reference@meta.data[, i])
  }
  reference@meta.data <- dplyr::select(reference@meta.data, -c("name"))


  for (i in seq(from = start.resolution, to = end.resolution, by = step)) {
    query <- Seurat::FindClusters(query, resolution = i)
    query_cluster_cells <- split(query@meta.data, query@meta.data$seurat_clusters)
    num.cells <- c()
    min.cells <- c()
    for (j in 1:length(query_cluster_cells)) {
      num.cells[j] <- nrow(query_cluster_cells[[j]])
    }
    min.cells <- min(num.cells)
    if(min.cells < min.cluster.cells)
    {
      query@meta.data <- query@meta.data[, -ncol(query@meta.data)]
      query@meta.data$seurat_clusters <- query@meta.data[, ncol(query@meta.data)]
      break
    }
  }
  query@meta.data <- dplyr::select(query@meta.data, -c("seurat_clusters"))
  query@meta.data$name <- "D2_"
  for (i in 4:(ncol(query@meta.data)-1)) {
    query@meta.data[, i] <- stringr::str_c(query$name,query@meta.data[, i])
  }
  query@meta.data <- dplyr::select(query@meta.data, -c("name"))

  #CCA embedding
  CCA <- Seurat::RunCCA(object1 = reference, object2 = query,num.cc = dim.cca)
  L2CCA<-Seurat::L2Dim(CCA,reduction = "cca")
  embedding<-L2CCA@reductions[["cca.l2"]]@cell.embeddings
  embedding<-data.frame(cbind(cells=rownames(embedding),embedding))

  #CCA cluster representation
  #reference
  reference_embedding_list <- list()
  for (x1 in 4:ncol(reference@meta.data)) {
    reference_cluster <- split(reference@meta.data, reference@meta.data[, x1])
    reference_means <- data.frame()
    for (i in 1:length(reference_cluster)) {
      cluster_name <- names(reference_cluster)[i]
      cluster_cell <- data.frame(cells = rownames(reference_cluster[[i]]))
      cluster_embedding <- dplyr::left_join(cluster_cell, embedding,by="cells")
      rownames(cluster_embedding)<-cluster_embedding[, 1]
      cluster_embedding <- cluster_embedding[, -1]
      cluster_embedding <- as.data.frame(lapply(cluster_embedding, as.numeric))
      cluster_vec <- colMeans(cluster_embedding)
      cluster_vec <- data.frame(cluster_vec)
      colnames(cluster_vec)[1] <- cluster_name
      if (i==1) {
        reference_means <- cluster_vec

      }
      else
        reference_means <- cbind(reference_means,cluster_vec)
    }
    reference_embedding_list[[x1-3]] <- reference_means
  }

  query_embedding_list <- list()
  for (x1 in 4:ncol(query@meta.data)) {
    query_cluster <- split(query@meta.data, query@meta.data[, x1])
    query_means <- data.frame()
    for (i in 1:length(query_cluster)) {
      cluster_name <- names(query_cluster)[i]
      cluster_cell <- data.frame(cells=rownames(query_cluster[[i]]))
      cluster_embedding <- dplyr::left_join(cluster_cell, embedding,by="cells")
      rownames(cluster_embedding) <- cluster_embedding[, 1]
      cluster_embedding <- cluster_embedding[, -1]
      cluster_embedding <- as.data.frame(lapply(cluster_embedding, as.numeric))
      cluster_vec <- colMeans(cluster_embedding)
      cluster_vec <- data.frame(cluster_vec)
      colnames(cluster_vec)[1] <- cluster_name
      if (i==1) {
        query_means <- cluster_vec

      }
      else
        query_means <- cbind(query_means, cluster_vec)
    }
    query_embedding_list[[x1-3]] <- query_means
  }


  #correlation of CCA cluster representation
  CCA_PCC_list <- list()
  for (x1 in 1:length(reference_embedding_list)) {
    PCC_list <- list()
    for (y1 in 1:length(query_embedding_list)) {
      CCA_PCC <- cor(query_embedding_list[[y1]], reference_embedding_list[[x1]])
      PCC_list[[y1]] <- CCA_PCC
    }
    CCA_PCC_list[[x1]] <- PCC_list
  }

  #Epsilon
  reference_all_means <- embedding[1:ncol(reference_matrix), ]
  reference_all_means <- reference_all_means[, -1]
  reference_all_means <- as.data.frame(lapply(reference_all_means, as.numeric))
  reference_all_means_vec <- colMeans(reference_all_means)
  query_all_means <- embedding[(ncol(reference_matrix)+1):nrow(embedding), ]
  query_all_means <- query_all_means[, -1]
  query_all_means <- as.data.frame(lapply(query_all_means, as.numeric))
  query_all_means_vec <- colMeans(query_all_means)
  Epsilon <- cor(reference_all_means_vec, query_all_means_vec)

  #match
  strongest_match_avg <- data.frame()
  for (x1 in 1:length(CCA_PCC_list)) {
    for (y1 in 1:length(CCA_PCC_list[[x1]])) {
      match_matrix <- Match(data = CCA_PCC_list[[x1]][[y1]], top = 1, threshold = Epsilon)
      strongest_match_PCC <- CCA_PCC_list[[x1]][[y1]][match_matrix == 2]
      strongest_match_avg[x1,y1] <- mean(strongest_match_PCC)
    }
  }

  #output
  name_row <- data.frame(res=seq(from = start.resolution, to = (ncol(reference@meta.data)-3)*step, by = step))
  name_row$name <- "D1_"
  name_row$res <- stringr::str_c(name_row$name, name_row$res)
  rownames(strongest_match_avg) <- name_row$res

  name_col <- data.frame(res=seq(from = start.resolution, to = (ncol(query@meta.data)-3)*step, by = step))
  name_col$name <- "D2_"
  name_col$res <- stringr::str_c(name_col$name, name_col$res)
  colnames(strongest_match_avg) <- name_col$res

  a <- as.matrix(strongest_match_avg)
  b <- which(a==max(a),arr.ind=T)
  b <- as.data.frame(b)
  c <- nrow(b)
  D1_res <- start.resolution + step*(b[c, 1]-1)
  D2_res <- start.resolution + step*(b[c, 2]-1)
  ClusterMatch_resolution <- list(Resolution_cor = strongest_match_avg, D1_ref_res = D1_res, D2_que_res =D2_res)

  return(ClusterMatch_resolution)
}

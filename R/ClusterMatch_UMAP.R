#' UMAP visualization
#' @param embedding cell embeddings
#' @param celltype a data.frame, the first column is the cell and the rest are the cell's label
#' @return a data.frame, the first column is the cell, the second and third columns are UMAP1 and UMAP2, and the rest are the cell's label
#' @importFrom umap umap.defaults umap
#' @export

ClusterMatch_UMAP <- function(embedding, celltype){
  custom.config <- umap::umap.defaults
  custom.config$random_state <- 123
  r_q.umap <- umap::umap(embedding)
  r_q <- data.frame(cell=rownames(r_q.umap$layout),r_q.umap$layout)
  r_q <- dplyr::left_join(r_q,celltype)
  return(r_q)
}






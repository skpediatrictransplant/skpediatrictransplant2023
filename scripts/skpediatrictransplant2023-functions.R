library(MASS)
library(stats)
library(lsa)
library(vegan)
library(FlowSOM)
library(variancePartition)
library(BiocParallel)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)


pca_function = function(data, features){
     pca.res = stats::prcomp(data %>% ungroup() %>% dplyr::select(all_of(features)), scale = T, center = T)
     pca.df= pca.res$x %>% data.frame() %>%
          bind_cols(data)
     return(pca.df)
}




lda_function = function(data,x,y){
     lda.res = MASS::lda(x, y,)
     lda.df = hsslda::makeAxes(data, co = lda.res$scaling)
     return(list(lda.df=lda.df, lda.res = lda.res))
}



umap_function = function(data, features){
     umap.res = uwot::umap(data %>% ungroup() %>% dplyr::select(all_of(features)), verbose= T)
     umap.df= umap.res %>% data.frame() %>%
          dplyr::rename(UMAP1 = X1, UMAP2= X2) %>%
          bind_cols(data)
     return(umap.df)
}


## function to calculate cosine similarity between each cell , within each sample
calculateCosineSimilarityOfEachSample = function(CD4_cells_sampleSplits){
     lapply(CD4_cells_sampleSplits, function(sampleSplit){
          if (nrow(sampleSplit) < 500){
               sampleSplit = sampleSplit  %>%
                    sample_n(nrow(sampleSplit))
          } else {
               sampleSplit =sampleSplit %>%
                    sample_n(500)
          }
          nm = unique(sampleSplit$sample) ## unique sample variable
          message('sample #: ', unique(sampleSplit$sample))
          cosine.res = cosine(sampleSplit %>% dplyr::select(all_of(features)) %>% t() %>% as.matrix())
          cosine.res[lower.tri(cosine.res, diag = T)] <- NA ## remove diagonal and lower half
          cosine.res = data.frame(cosine.res)
          colnames(cosine.res) = sampleSplit$cell.id
          rownames(cosine.res) = sampleSplit$cell.id
          
          cosine.res_pairwise = cosine.res %>%
               data.frame() %>%
               rownames_to_column('cell.id2') %>%
               gather(key = 'cell.id', value = 'cosine', -cell.id2) %>%
               na.omit() %>%
               left_join(CD4_cells %>% dplyr::select(cell.id, celltype), by = 'cell.id') %>% 
               left_join(CD4_cells %>% dplyr::select(cell.id2 = cell.id, celltype2 = celltype), by = 'cell.id2') %>%
               group_by(celltype,celltype2) %>%
               summarize(cosine = mean(cosine),) %>%
               mutate(sample = nm)
          return(cosine.res_pairwise)
     })
}




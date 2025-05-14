
wind <- pexp$specRegion.inds

specRegion = pexp$specRegion
ppmRegion = pexp$ppmRegion
n_components_to_remove <- 1

remove_background_pca <- function(specRegion, n_components_to_remove = 1) {
  # Center the data
  centered <- scale(specRegion, center = TRUE, scale = FALSE)
  
  # Perform PCA
  pca <- prcomp(centered, center = F, scale. = F)
  
  # Reconstruct the background using leading PCs (e.g., PC1)
  background_estimate <- pca$x[, 1:n_components_to_remove] %*% 
                         t(pca$rotation[, 1:n_components_to_remove])
  
  # Subtract background from centered data
  corrected <- centered - background_estimate
  return(corrected)
  # Optionally: add back the original mean to maintain absolute intensities
  # corrected <- lapply(1:nrow(corrected), function(x){
  #   corrected[x, ] + attr(centered, "scaled:center")
  # }) %>% do.call(rbind, .)
  # 
  
}

corrected <- remove_background_pca(specRegion, n_components_to_remove=1)
simplePlot(corrected, ppmRegion)
simplePlot(specRegion, ppmRegion)
simplePlot(abs(background_estimate), ppmRegion)

# This source file includes separate functions for the following purposes
# 1) Standardize pixel intensity
# 2) get indices for convoluion
# 3) convolution with kernel or function
# 4) get image determinant of Hessian
# 5) get pixel values for image
# 6) get image local maxima


############################### Standardize pixel intensity ##############################

standardize <- function(
  image_obj, 
  ...)
  
{
  img <- image_obj
  
  # handle missing values
  if (any(is.na(img))) {
    warning(paste0("The image has ", ncell(img), " missing values (", 
                   sum(is.na(img))/ncell(img), "%)."))
  }
  
  i_min <- min(img, ...) # get minimum intensity
  i_max <- max(img, ...) # get maximum intensity
  i_range <- i_max - i_min # get intensity range
  
  return((img - i_min)/i_range) # return standardized image
}



############################### get indices for convoluion ###############################

conv_indices <- function(
  image_obj, 
  kernel_window)
  
{
  
  wind <- kernel_window 
  mrg <- wind %/% 2 # use the quotent of wind/2 for use as a margin and sd of gaussian PSF 
  img <- image_obj
  
  #### x axis ####
  
  x_coord <- 1:nrow(img) # get initial coordinates along horizontal axis
  rem <- length(x_coord) %% wind # get remainder of (x_coord length)/wind for margin cliping
  if (rem > 0) { 
    x_coord <- tail(x_coord, n = -rem) # remove margin of size 'rem' from the right
  }
  x_coord <- tail(x_coord, n = -mrg) 
  x_coord <- head(x_coord, n = -mrg)
  
  #### y axis ####
  
  y_coord <- 1:ncol(img)
  rem <- length(y_coord) %% wind
  if (rem > 0) {
    y_coord <- tail(y_coord, n = -rem)
  }
  
  y_coord <- tail(y_coord, n = -mrg)
  y_coord <- head(y_coord, n = -mrg)
  
  #### make data frame with indices ####
  return(expand.grid(x = x_coord, y = y_coord))
}  



########################### convolution with kernel or function ##########################


conv_image <- function(image_obj = crop1, 
                       kernel_window. = 3,
                       kernel_matrix = kernel_matrix) 
{
  img <- image_obj
  mrg <- kernel_window. %/% 2 # get kernel margins
  ind <- conv_indices(image_obj = img, kernel_window = kernel_window.) # get conv indices
  
  for (i in 1:nrow(ind)) {
    x <- ind[i, 1, drop = T]
    y <- ind[i, 2, drop = T]
    original_kernel <- image_obj[(x - mrg):(x + mrg), (y - mrg):(y + mrg)]
    
    img[x, y] <- sum(original_kernel * kernel_matrix)
  }
  return(img)
}

########################### get image determinant of Hessian #############################


Hess_det <- function(image_obj) {
  if (!require(imager)) {
    message("Installing the following required package: imager")
    install.packages("imager")
    library(imager)
  }
  
  hess <- imhessian(as.cimg(image_obj))
  det_hess <- with(hess, (xx*yy - xy^2))
  return(as.matrix(det_hess))
}

############################## get pixel values for image ################################


pixel_values <- function(image_obj, xy_coord) {
  
  val <- vector(length = nrow(xy_coord))
  
  for (i in 1:nrow(xy_coord)) {
    x <- xy_coord[i, 1, drop = T]
    y <- xy_coord[i, 2, drop = T]
    val[i] <- image_obj[x, y]
  }
  
  return(val)
}

################################# get image local maxima sumary #################################


get_local_max <- function(image_obj, 
                          kernel_window,
                          thresh_coef = 1.3) {
  
  img <- image_obj
  det_hess <- Hess_det(img)
  
  wind <- kernel_window
  mrg <- wind %/% 2
  
  
  indices <- conv_indices(img, wind)
  intensity <- vector(length = nrow(indices))
  kernel_mean_intensity <- vector(length = nrow(indices))
  is_max <- vector(length = nrow(indices))
  DHess <- vector(length = nrow(indices))
  mean_x_DHess <- vector(length = nrow(indices))
  
  table. <- data.frame(indices, 
                       intensity, 
                       kernel_mean_intensity,
                       is_max, 
                       DHess,
                       mean_x_DHess)

  for (i in 1:nrow(table.)) {
    x <- table.[i, 1, drop = T]
    y <- table.[i, 2, drop = T]
    
    img_kernel <- img[(x - mrg):(x + mrg), 
                      (y - mrg):(y + mrg)]
    kernel_mean <- mean(img_kernel)
    
    det_hess. <- det_hess[x, y]
    
    table.[i, "intensity"] <- img[x, y]
    
    table.[i, "kernel_mean_intensity"] <- kernel_mean
    
    table.[i, "is_max"] <- img[x, y] >= max(img_kernel)
    
    table.[i, "DHess"] <- det_hess.
    
    table.[i, "mean_x_DHess"] <- kernel_mean *  det_hess. 
  }
  

  background_mean <- mean(table.[!table.[, "is_max", drop = T], "intensity", drop = T])
  background_sd <- sd(table.[!table.[, "is_max", drop = T], "intensity", drop = T])
  
  significance <- thresh_coef * (table.[, "intensity", drop = T] - background_mean)/background_sd
  is_sig <- table.[, "intensity", drop = T] > significance 
  is_max_sig <- table.[, "is_max", drop = T] & is_sig
  
  table. <- cbind(table., 
                  significance, 
                  is_sig, 
                  is_max_sig)
  
  local_maxima_only_table <- table.[table.[, "is_max", drop = T], -c(5, 9)]
  
  tables <- list(length = 2)
  
  tables[[1]] <- table.
  tables[[2]] <- local_maxima_only_table
  
  names(tables) <- c("image_summary", "image_summary_LMonly")
  
  return(tables)
}
  


######################### identify particles from local maxima ###########################
################################# make image plot #################################
image_obj = crop1_bx
point_indices <- particle_indices 
  

print_image <- function(image_obj,
                        add_label = F,
                        point_indices = NULL,
                        rsn = 150,
                        file_name,
                        quant,
                        ...) 
{
  img <- as.cimg(image_obj)
  fn <- paste0(file_name, ".tif")
  
  tiff(fn, res = rsn)
  par(mar = c(0, 0, 0, 0), ...)
  plot(img)
  
  if (!add_label) {
    dev.off() 
  } else {
    Dhess_mean <- image_summary_LMonly$mean_x_DHess
    cutoff <- quantile(Dhess_mean, quant)
    particle_indices <- (image_summary_LMonly[Dhess_mean > cutoff, 1:2])
    
    tiff(fn, res = rsn)
    par(mar = c(0, 0, 0, 0))
    plot(img)
    with(particle_indices, points(x, y, col = "blue", lwd = 1))
    dev.off()
  }
  
}




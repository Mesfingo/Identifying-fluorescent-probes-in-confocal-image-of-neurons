# do particle detection

# 1) get image data, change to 2D matrix and standardize between 0 and 1 
library("raster")
library("rgdal")
library("imager")

image_org <- raster("c1_green_C0_T4.ome.tif")
crop1 <- (image_org[600:725, 450:575, drop = F])
crop1 <- as.matrix(crop1)
crop1 <- standardize(crop1)

# 2) Low pass (box average) filter
kernel_window <- 3
kernel_coef <- 1/(kernel_window)^2
kernel_matrix <- matrix(rep(kernel_coef, (kernel_window)^2),
                        nrow = kernel_window,
                        ncol = kernel_window)

crop1_bx <- conv_image(image_obj = crop1, kernel_matrix = kernel_matrix)

# 3) hi pass (Gaussian) filter

Gau_kernl <- matrix(c( 0.0585, 0.0965, 0.0585, 
                       0.0965, 0.1592, 0.0965,
                       0.0585, 0.0965, 0.0585),
                    nrow = 3, byrow = T)

crop1_bx_Gu <- conv_image(image_obj = crop1_bx, kernel_matrix = Gau_kernl)

# 4) Laplacial of Gaussian filter

LoG_kernl <- matrix(c( 0.0000, -0.0965,  0.0000, 
                       -0.0965, -0.3183, -0.0965,
                       0.0000, -0.0965,  0.0000),
                    nrow = 3, byrow = T)

crop1_bx_Gu_LoG <- conv_image(image_obj = crop1_bx_Gu, kernel_matrix = LoG_kernl)

crop1_bx_Gu_LoG <- crop1_bx_Gu_LoG + crop1

# 5) Detect local maxima on filteres image

detection_summary <- get_local_max(crop1_bx_Gu_LoG, kernel_window = 3)

image_summary <- detection_summary[[1]]
image_summary_LMonly <- detection_summary[[2]]



#unloadNamespace("rasterVis")
#unloadNamespace("raster")
Dhess_mean <- image_summary_LMonly$mean_x_DHess
hist(Dhess_mean, breaks = 20)
hist(Dhess_mean, breaks = 25)
hist(Dhess_mean, breaks = 20)
hist(Dhess_mean, breaks = 20)

#par(col = "#8b0000", , cex = 1)
#boxplot(Dhess_mean, outcol = "#8b0000", pch = 21, bty = "l")

#xlim=c(0,150), ylim=c(0,.05)
hist(Dhess_mean, breaks = 40, 
     freq = F, prob = T)
lines(density(Dhess_mean, bw = 1), col = "#8b0000")

hist(Dhess_mean, prob = TRUE, breaks = 100, col = "#8b0000")            # prob=TRUE for probabilities not counts
lines(density(Dhess_mean))             # add a density estimate with defaults
lines(density(Dhess_mean, adjust = 3), col = "red")


#plot(ecdf(Dhess_mean))
#h <- hist(Dhess_mean, breaks = 50, col = "red", prob = TRUE)
#h$counts <- cumsum(h$counts)
#plot(h)

png("Spot Classification Response3.png", res = 100,
    height = 600,
    width = 800)
plot(sort(Dhess_mean), 
     1:length(Dhess_mean),
     xlab = "Spot Classification Response",
     ylab = "Local Maxima Spot Count",
     col = "#8b0000",
     bty = "l",
     pch = 20)
abline(v = 0.01, lty = 2, col = "#ff4500")
abline(v = quantile(Dhess_mean, 0.99), 
      lty = 4, col = "#00008b")
dev.off()

cutoff <- quantile(Dhess_mean, 0.99)

# get indices of particles above treshold

particle_indices <- (image_summary_LMonly[Dhess_mean > cutoff, 1:2])


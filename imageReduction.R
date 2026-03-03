# For image reading and manipulation
library(jpeg)        # for reading JPG images
library(grid)        # for grid.raster plotting

img <- readJPEG("i.jpg")  # reads the image into a 3D array (rows x cols x channels)
dim(img)                  # check dimensions: height x width x 3 (RGB)

# Get image dimensions

# We need to convert the image into a 2D matrix where each row is a pixel and each column is a color channel (R, G, B)

img_dim <- dim(img)
height <- img_dim[1]
width  <- img_dim[2]

# Flatten the image to a 2D matrix (pixels x RGB)
img_matrix <- matrix(
  c(img[,,1], img[,,2], img[,,3]),
  ncol = 3
)

#Each row in img_matrix now represents the RGB values of one pixel.

# total pixels : height x width 
# In most images, each channel (R, G, B) can be represented in 8-bit integer form: 0 to 255.
#However jpeg scales 0-1 , original value/255. So e.g pure red : [225,0,0]

set.seed(123)  # for reproducibility

kmeans_results <- list()
cluster_sizes <- c(100)

#we want to group pixels into 2 groups first.

for (k in cluster_sizes) {
  kmeans_results[[as.character(k)]] <- kmeans(img_matrix, centers = k, iter.max = 20)
}

#so it calculate a 3d euclidean distance x,y,z 
#the assign the point to the each cluster based on closeness to a 3d centroid

compressed_images <- list()

for (k in cluster_sizes) {
  km <- kmeans_results[[as.character(k)]]
  
  # Get cluster centers
  centers <- km$centers
  
  # Replace each pixel with its centroid
  clustered_pixels <- centers[km$cluster, ]
  
  # Reshape back to image dimensions
  compressed_img <- array(clustered_pixels, dim = c(height, width, 3))
  
  compressed_images[[as.character(k)]] <- compressed_img
}

# Set up a 1-row layout for all cluster images
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, length(cluster_sizes))))

for (i in seq_along(cluster_sizes)) {
  k <- cluster_sizes[i]
  img_to_plot <- compressed_images[[as.character(k)]]
  
  # Position in layout
  print(
    grid.raster(img_to_plot),
    vp = viewport(layout.pos.row = 1, layout.pos.col = i)
  )
}




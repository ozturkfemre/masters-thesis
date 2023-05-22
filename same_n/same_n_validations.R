library(mvtnorm)
library(tidyverse)
library(cluster)
library(NbClust)
library(fpc)
library(factoextra)
library(fossil)
####################


generateGaussianData <- function(n, center, sigma, label) {
  data = rmvnorm(n, mean = center, sigma = sigma)
  data = data.frame(data)
  data = data %>% mutate(class=factor(label))
  data
}


################################################################################
# //////////////////////////////////////////////////////////////////////////// #
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ###
#### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ####
################################################################################

#############################################
### Validations for 2 clustered data sets ###
#############################################


##############
### Case A ###
##############


### A1 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 50, 50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 0, 0 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # all data
  data <- bind_rows(data1, data2)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)

# Since all of the methdos offered 2 clusters, the same results will be generated. I will do this part only once.


km_data <- eclust(df, "kmeans", k = 2, nstart = 25, graph = F)


rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))


cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi

two_k_same_n_valid <- list()

two_k_same_n_valid$silhouette <- data.frame()
two_k_same_n_valid$ch <- data.frame()
two_k_same_n_valid$db <- data.frame()
two_k_same_n_valid$ascm <- data.frame()
 
two_k_same_n_valid$silhouette[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$silhouette[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi

two_k_same_n_valid$ch[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$ch[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi

two_k_same_n_valid$db[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$db[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi

two_k_same_n_valid$ascm[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$ascm[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi


names(two_k_same_n_valid$silhouette) <- c("rand", "mvi")
names(two_k_same_n_valid$ch) <- c("rand", "mvi")
names(two_k_same_n_valid$db) <- c("rand", "mvi")
names(two_k_same_n_valid$ascm) <- c("rand", "mvi")

### A2 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 50, 50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 0, 0 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # all data
  data <- bind_rows(data1, data2)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)

# Since all of the methdos offered 2 clusters, the same results will be generated. I will do this part only once.


km_data <- eclust(df, "kmeans", k = 2, nstart = 25, graph = F)



two_k_same_n_valid$silhouette[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$silhouette[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi

two_k_same_n_valid$ch[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$ch[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi

two_k_same_n_valid$db[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$db[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi

two_k_same_n_valid$ascm[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$ascm[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi



### A3 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 50, 50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 0, 0 )
  sigma <- matrix(c(75,0,0,75), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # all data
  data <- bind_rows(data1, data2)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)

# Since all of the methdos offered 2 clusters, the same results will be generated. I will do this part only once.


km_data <- eclust(df, "kmeans", k = 2, nstart = 25, graph = F)



two_k_same_n_valid$silhouette[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$silhouette[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi

two_k_same_n_valid$ch[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$ch[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi

two_k_same_n_valid$db[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$db[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi

two_k_same_n_valid$ascm[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$ascm[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi


##############
### Case B ###
##############

### B1 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 50, 50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 25, 25 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # all data
  data <- bind_rows(data1, data2)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)

# Since all of the methdos offered 2 clusters, the same results will be generated. I will do this part only once.


km_data <- eclust(df, "kmeans", k = 2, nstart = 25, graph = F)



two_k_same_n_valid$silhouette[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$silhouette[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi

two_k_same_n_valid$ch[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$ch[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi

two_k_same_n_valid$db[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$db[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi

two_k_same_n_valid$ascm[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$ascm[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi



### B2 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 50, 50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 25, 25 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # all data
  data <- bind_rows(data1, data2)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)

# Since all of the methdos offered 2 clusters, the same results will be generated. I will do this part only once.


km_data <- eclust(df, "kmeans", k = 2, nstart = 25, graph = F)



two_k_same_n_valid$silhouette[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$silhouette[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi

two_k_same_n_valid$ch[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$ch[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi

two_k_same_n_valid$db[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$db[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi

two_k_same_n_valid$ascm[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$ascm[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi


### B3 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 50, 50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 25, 25 )
  sigma <- matrix(c(75,0,0,75), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # all data
  data <- bind_rows(data1, data2)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)

# Since all of the methdos offered 2 clusters, the same results will be generated. I will do this part only once.


km_data <- eclust(df, "kmeans", k = 2, nstart = 25, graph = F)



two_k_same_n_valid$silhouette[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$silhouette[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi

two_k_same_n_valid$ch[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$ch[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi

two_k_same_n_valid$db[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$db[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi

two_k_same_n_valid$ascm[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(km_data$cluster))
two_k_same_n_valid$ascm[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), km_data$cluster)$vi



################################################################################
# //////////////////////////////////////////////////////////////////////////// #
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ###
#### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ####
################################################################################

#############################################
### Validations for 2 clustered data sets ###
#############################################


##############
### Case A ###
##############


### A1 ###
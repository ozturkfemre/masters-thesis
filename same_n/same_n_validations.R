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

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 50,100 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 25,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( 50,0 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # all data
  data <- bind_rows(data1, data2, data3)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)


three_k_same_n_valid <- list()
three_k_same_n_valid$silhouette <- data.frame()
three_k_same_n_valid$ch <- data.frame()
three_k_same_n_valid$db <- data.frame()
three_k_same_n_valid$ascm <- data.frame()


kmsil_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)



three_k_same_n_valid$silhouette[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
three_k_same_n_valid$silhouette[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

three_k_same_n_valid$ch[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
three_k_same_n_valid$ch[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

three_k_same_n_valid$db[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
three_k_same_n_valid$db[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

three_k_same_n_valid$ascm[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
three_k_same_n_valid$ascm[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi


names(three_k_same_n_valid$silhouette) <- c("rand", "mvi")
names(three_k_same_n_valid$ch) <- c("rand", "mvi")
names(three_k_same_n_valid$db) <- c("rand", "mvi")
names(three_k_same_n_valid$ascm) <- c("rand", "mvi")



### A2 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 50,100 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 25,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( 50,0 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # all data
  data <- bind_rows(data1, data2, data3)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)



kmsil_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)



three_k_same_n_valid$silhouette[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
three_k_same_n_valid$silhouette[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

three_k_same_n_valid$ch[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
three_k_same_n_valid$ch[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

three_k_same_n_valid$db[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
three_k_same_n_valid$db[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

three_k_same_n_valid$ascm[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
three_k_same_n_valid$ascm[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi


### A3 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 50,100 )
  sigma <- matrix(c(100,0,0,100), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 25,50 )
  sigma <- matrix(c(80,0,0,80), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( 50,0 )
  sigma <- matrix(c(60,0,0,60), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # all data
  data <- bind_rows(data1, data2, data3)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)



kmsil_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)



three_k_same_n_valid$silhouette[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
three_k_same_n_valid$silhouette[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

three_k_same_n_valid$ch[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
three_k_same_n_valid$ch[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

three_k_same_n_valid$db[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
three_k_same_n_valid$db[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

three_k_same_n_valid$ascm[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
three_k_same_n_valid$ascm[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi



##############
### Case B ###
##############

### B1 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 75,75 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 50,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( 25,75 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # all data
  data <- bind_rows(data1, data2, data3)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)


kmsil_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)



three_k_same_n_valid$silhouette[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
three_k_same_n_valid$silhouette[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

three_k_same_n_valid$ch[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
three_k_same_n_valid$ch[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

three_k_same_n_valid$db[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
three_k_same_n_valid$db[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

three_k_same_n_valid$ascm[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
three_k_same_n_valid$ascm[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi


### B2 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 75,75 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 50,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( 25,75 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # all data
  data <- bind_rows(data1, data2, data3)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)


kmsil_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)



three_k_same_n_valid$silhouette[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
three_k_same_n_valid$silhouette[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

three_k_same_n_valid$ch[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
three_k_same_n_valid$ch[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

three_k_same_n_valid$db[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
three_k_same_n_valid$db[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

three_k_same_n_valid$ascm[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
three_k_same_n_valid$ascm[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi



### B3 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 75,75 )
  sigma <- matrix(c(100,0,0,100), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 50,50 )
  sigma <- matrix(c(80,0,0,80), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( 25,75 )
  sigma <- matrix(c(60,0,0,60), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # all data
  data <- bind_rows(data1, data2, data3)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)


kmsil_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)



three_k_same_n_valid$silhouette[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
three_k_same_n_valid$silhouette[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

three_k_same_n_valid$ch[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
three_k_same_n_valid$ch[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

three_k_same_n_valid$db[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
three_k_same_n_valid$db[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

three_k_same_n_valid$ascm[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
three_k_same_n_valid$ascm[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi



################################################################################
# //////////////////////////////////////////////////////////////////////////// #
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ###
#### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ####
################################################################################

#############################################
### Validations for 4 clustered data sets ###
#############################################


##############
### Case A ###
##############


### A1 ###


dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,25 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 0,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( 50,75 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # all data
  data <- bind_rows(data1, data2, data3, data4)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)



four_k_same_n_valid <- list()
four_k_same_n_valid$silhouette <- data.frame()
four_k_same_n_valid$ch <- data.frame()
four_k_same_n_valid$db <- data.frame()
four_k_same_n_valid$ascm <- data.frame()


kmsil_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 5, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)



four_k_same_n_valid$silhouette[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
four_k_same_n_valid$silhouette[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

four_k_same_n_valid$ch[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
four_k_same_n_valid$ch[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

four_k_same_n_valid$db[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
four_k_same_n_valid$db[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

four_k_same_n_valid$ascm[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
four_k_same_n_valid$ascm[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi


names(four_k_same_n_valid$silhouette) <- c("rand", "mvi")
names(four_k_same_n_valid$ch) <- c("rand", "mvi")
names(four_k_same_n_valid$db) <- c("rand", "mvi")
names(four_k_same_n_valid$ascm) <- c("rand", "mvi")



### A2 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,25 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 0,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( 50,75 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # all data
  data <- bind_rows(data1, data2, data3, data4)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)



kmsil_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 5, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)



four_k_same_n_valid$silhouette[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
four_k_same_n_valid$silhouette[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

four_k_same_n_valid$ch[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
four_k_same_n_valid$ch[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

four_k_same_n_valid$db[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
four_k_same_n_valid$db[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

four_k_same_n_valid$ascm[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
four_k_same_n_valid$ascm[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi



### A3 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,25 )
  sigma <- matrix(c(100,0,0,100), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 0,50 )
  sigma <- matrix(c(80,0,0,80), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( 50,75 )
  sigma <- matrix(c(60,0,0,60), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(40,0,0,40), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # all data
  data <- bind_rows(data1, data2, data3, data4)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)



kmsil_data <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 5, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)



four_k_same_n_valid$silhouette[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
four_k_same_n_valid$silhouette[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

four_k_same_n_valid$ch[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
four_k_same_n_valid$ch[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

four_k_same_n_valid$db[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
four_k_same_n_valid$db[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

four_k_same_n_valid$ascm[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
four_k_same_n_valid$ascm[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi



##############
### Case B ###
##############

### B1 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 100,75 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 75,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( 50,75 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 25,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # all data
  data <- bind_rows(data1, data2, data3, data4)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)



kmsil_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)



four_k_same_n_valid$silhouette[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
four_k_same_n_valid$silhouette[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

four_k_same_n_valid$ch[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
four_k_same_n_valid$ch[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

four_k_same_n_valid$db[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
four_k_same_n_valid$db[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

four_k_same_n_valid$ascm[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
four_k_same_n_valid$ascm[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi


### B2 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 100,75 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 75,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( 50,75 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 25,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # all data
  data <- bind_rows(data1, data2, data3, data4)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)


kmsil_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)



four_k_same_n_valid$silhouette[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
four_k_same_n_valid$silhouette[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

four_k_same_n_valid$ch[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
four_k_same_n_valid$ch[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

four_k_same_n_valid$db[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
four_k_same_n_valid$db[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

four_k_same_n_valid$ascm[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
four_k_same_n_valid$ascm[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi



### B3 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 100,75 )
  sigma <- matrix(c(100,0,0,100), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 75,50 )
  sigma <- matrix(c(80,0,0,80), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( 50,75 )
  sigma <- matrix(c(60,0,0,60), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 25,50 )
  sigma <- matrix(c(40,0,0,40), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # all data
  data <- bind_rows(data1, data2, data3, data4)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)


kmsil_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)



four_k_same_n_valid$silhouette[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
four_k_same_n_valid$silhouette[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

four_k_same_n_valid$ch[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
four_k_same_n_valid$ch[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

four_k_same_n_valid$db[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
four_k_same_n_valid$db[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

four_k_same_n_valid$ascm[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
four_k_same_n_valid$ascm[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi




################################################################################
# //////////////////////////////////////////////////////////////////////////// #
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ###
#### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ####
################################################################################

#############################################
### Validations for 5 clustered data sets ###
#############################################


##############
### Case A ###
##############


### A1 ###


dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 200,25 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 150,100 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( 50,100 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 50,25 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( -50,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)



five_k_same_n_valid <- list()
five_k_same_n_valid$silhouette <- data.frame()
five_k_same_n_valid$ch <- data.frame()
five_k_same_n_valid$db <- data.frame()
five_k_same_n_valid$ascm <- data.frame()


kmsil_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 5, nstart = 25, graph = F)



five_k_same_n_valid$silhouette[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
five_k_same_n_valid$silhouette[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

five_k_same_n_valid$ch[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
five_k_same_n_valid$ch[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

five_k_same_n_valid$db[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
five_k_same_n_valid$db[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

five_k_same_n_valid$ascm[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
five_k_same_n_valid$ascm[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi


names(five_k_same_n_valid$silhouette) <- c("rand", "mvi")
names(five_k_same_n_valid$ch) <- c("rand", "mvi")
names(five_k_same_n_valid$db) <- c("rand", "mvi")
names(five_k_same_n_valid$ascm) <- c("rand", "mvi")



### A2 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 200,25 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 150,100 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( 50,100 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 50,25 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( -50,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)


kmsil_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 5, nstart = 25, graph = F)



five_k_same_n_valid$silhouette[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
five_k_same_n_valid$silhouette[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

five_k_same_n_valid$ch[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
five_k_same_n_valid$ch[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

five_k_same_n_valid$db[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
five_k_same_n_valid$db[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

five_k_same_n_valid$ascm[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
five_k_same_n_valid$ascm[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi




### A3 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 200,25 )
  sigma <- matrix(c(100,0,0,100), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 150,100 )
  sigma <- matrix(c(80,0,0,80), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( 50,100 )
  sigma <- matrix(c(60,0,0,60), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 50,25 )
  sigma <- matrix( c(40,0,0,40), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( -50,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)



kmsil_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 5, nstart = 25, graph = F)



five_k_same_n_valid$silhouette[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
five_k_same_n_valid$silhouette[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

five_k_same_n_valid$ch[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
five_k_same_n_valid$ch[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

five_k_same_n_valid$db[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
five_k_same_n_valid$db[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

five_k_same_n_valid$ascm[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
five_k_same_n_valid$ascm[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi




##############
### Case B ###
##############

### B1 ###


dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 150,25 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 100,100 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( 50,100 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 50,25 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( -50,100 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)


kmsil_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 5, nstart = 25, graph = F)



five_k_same_n_valid$silhouette[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
five_k_same_n_valid$silhouette[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

five_k_same_n_valid$ch[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
five_k_same_n_valid$ch[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

five_k_same_n_valid$db[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
five_k_same_n_valid$db[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

five_k_same_n_valid$ascm[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
five_k_same_n_valid$ascm[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi




### B2 ###


dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 150,25 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 100,100 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( 50,100 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 50,25 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( -50,100 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)



kmsil_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 5, nstart = 25, graph = F)



five_k_same_n_valid$silhouette[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
five_k_same_n_valid$silhouette[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

five_k_same_n_valid$ch[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
five_k_same_n_valid$ch[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

five_k_same_n_valid$db[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
five_k_same_n_valid$db[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

five_k_same_n_valid$ascm[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
five_k_same_n_valid$ascm[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi




### B3 ###


dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( 150,25 )
  sigma <- matrix(c(100,0,0,100), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( 100,100 )
  sigma <- matrix(c(80,0,0,80), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( 50,100 )
  sigma <- matrix(c(60,0,0,60), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 50,25 )
  sigma <- matrix( c(40,0,0,40), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( -50,100 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)



kmsil_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 4, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 5, nstart = 25, graph = F)



five_k_same_n_valid$silhouette[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
five_k_same_n_valid$silhouette[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

five_k_same_n_valid$ch[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
five_k_same_n_valid$ch[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

five_k_same_n_valid$db[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
five_k_same_n_valid$db[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

five_k_same_n_valid$ascm[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
five_k_same_n_valid$ascm[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi



################################################################################
# //////////////////////////////////////////////////////////////////////////// #
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ###
#### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ####
################################################################################

#############################################
### Validations for 5 clustered data sets ###
#############################################


##############
### Case A ###
##############


### A1 ###


dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,100 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( -150,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( -50,0 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 200,0 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # cluster 6
  n <- 50
  center <- c( 200,100 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data6 <- generateGaussianData(n, center, sigma, 6)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5, data6)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)


six_k_same_n_valid <- list()
six_k_same_n_valid$silhouette <- data.frame()
six_k_same_n_valid$ch <- data.frame()
six_k_same_n_valid$db <- data.frame()
six_k_same_n_valid$ascm <- data.frame()


kmsil_data <- eclust(df, "kmeans", k = 9, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 5, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)



six_k_same_n_valid$silhouette[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
six_k_same_n_valid$silhouette[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

six_k_same_n_valid$ch[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
six_k_same_n_valid$ch[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

six_k_same_n_valid$db[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
six_k_same_n_valid$db[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

six_k_same_n_valid$ascm[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
six_k_same_n_valid$ascm[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi


names(six_k_same_n_valid$silhouette) <- c("rand", "mvi")
names(six_k_same_n_valid$ch) <- c("rand", "mvi")
names(six_k_same_n_valid$db) <- c("rand", "mvi")
names(six_k_same_n_valid$ascm) <- c("rand", "mvi")



### A2 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,100 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( -150,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( -50,0 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 200,0 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # cluster 6
  n <- 50
  center <- c( 200,100 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data6 <- generateGaussianData(n, center, sigma, 6)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5, data6)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)


kmsil_data <- eclust(df, "kmeans", k = 9, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 9, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 7, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)



six_k_same_n_valid$silhouette[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
six_k_same_n_valid$silhouette[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

six_k_same_n_valid$ch[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
six_k_same_n_valid$ch[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

six_k_same_n_valid$db[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
six_k_same_n_valid$db[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

six_k_same_n_valid$ascm[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
six_k_same_n_valid$ascm[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi



### A3 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,100 )
  sigma <- matrix(c(100,0,0,100), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( -150,50 )
  sigma <- matrix(c(80,0,0,80), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( -50,0 )
  sigma <- matrix(c(60,0,0,60), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 200,0 )
  sigma <- matrix(c(40,0,0,40), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(70,0,0,70), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # cluster 6
  n <- 50
  center <- c( 200,100 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data6 <- generateGaussianData(n, center, sigma, 6)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5, data6)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)


kmsil_data <- eclust(df, "kmeans", k = 9, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 7, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)



six_k_same_n_valid$silhouette[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
six_k_same_n_valid$silhouette[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

six_k_same_n_valid$ch[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
six_k_same_n_valid$ch[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

six_k_same_n_valid$db[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
six_k_same_n_valid$db[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

six_k_same_n_valid$ascm[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
six_k_same_n_valid$ascm[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi



### B1 ###


dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,75 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( -150,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( -50,25 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 200,25 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # cluster 6
  n <- 50
  center <- c( 200,75 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data6 <- generateGaussianData(n, center, sigma, 6)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5, data6)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)

kmsil_data <- eclust(df, "kmeans", k = 9, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 7, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)



six_k_same_n_valid$silhouette[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
six_k_same_n_valid$silhouette[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

six_k_same_n_valid$ch[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
six_k_same_n_valid$ch[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

six_k_same_n_valid$db[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
six_k_same_n_valid$db[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

six_k_same_n_valid$ascm[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
six_k_same_n_valid$ascm[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi



### B2 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,75 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( -150,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( -50,25 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 200,25 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # cluster 6
  n <- 50
  center <- c( 200,75 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data6 <- generateGaussianData(n, center, sigma, 6)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5, data6)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)



kmsil_data <- eclust(df, "kmeans", k = 9, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 7, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)



six_k_same_n_valid$silhouette[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
six_k_same_n_valid$silhouette[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

six_k_same_n_valid$ch[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
six_k_same_n_valid$ch[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

six_k_same_n_valid$db[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
six_k_same_n_valid$db[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

six_k_same_n_valid$ascm[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
six_k_same_n_valid$ascm[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi



### B3 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,75 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( -150,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( -50,25 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 200,25 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # cluster 6
  n <- 50
  center <- c( 200,75 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data6 <- generateGaussianData(n, center, sigma, 6)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5, data6)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)



kmsil_data <- eclust(df, "kmeans", k = 9, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 7, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)



six_k_same_n_valid$silhouette[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
six_k_same_n_valid$silhouette[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

six_k_same_n_valid$ch[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
six_k_same_n_valid$ch[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

six_k_same_n_valid$db[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
six_k_same_n_valid$db[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

six_k_same_n_valid$ascm[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
six_k_same_n_valid$ascm[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi



################################################################################
# //////////////////////////////////////////////////////////////////////////// #
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ###
#### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ####
################################################################################

#############################################
### Validations for 7 clustered data sets ###
#############################################


##############
### Case A ###
##############



### A1 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,100 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( -150,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( -50,0 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 200,0 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # cluster 6
  n <- 50
  center <- c( 200,100 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data6 <- generateGaussianData(n, center, sigma, 6)
  # cluster 7
  n <- 50
  center <- c( -150,-50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data7 <- generateGaussianData(n, center, sigma, 7)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5, data6, data7)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)




seven_k_same_n_valid <- list()
seven_k_same_n_valid$silhouette <- data.frame()
seven_k_same_n_valid$ch <- data.frame()
seven_k_same_n_valid$db <- data.frame()
seven_k_same_n_valid$ascm <- data.frame()


kmsil_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 10, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 7, nstart = 25, graph = F)



seven_k_same_n_valid$silhouette[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
seven_k_same_n_valid$silhouette[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

seven_k_same_n_valid$ch[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
seven_k_same_n_valid$ch[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

seven_k_same_n_valid$db[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
seven_k_same_n_valid$db[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

seven_k_same_n_valid$ascm[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
seven_k_same_n_valid$ascm[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi


names(seven_k_same_n_valid$silhouette) <- c("rand", "mvi")
names(seven_k_same_n_valid$ch) <- c("rand", "mvi")
names(seven_k_same_n_valid$db) <- c("rand", "mvi")
names(seven_k_same_n_valid$ascm) <- c("rand", "mvi")



### A2 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,100 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( -150,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( -50,0 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 200,0 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # cluster 6
  n <- 50
  center <- c( 200,100 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data6 <- generateGaussianData(n, center, sigma, 6)
  # cluster 7
  n <- 50
  center <- c( -150,-50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data7 <- generateGaussianData(n, center, sigma, 7)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5, data6, data7)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)


kmsil_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 8, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 7, nstart = 25, graph = F)



seven_k_same_n_valid$silhouette[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
seven_k_same_n_valid$silhouette[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

seven_k_same_n_valid$ch[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
seven_k_same_n_valid$ch[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

seven_k_same_n_valid$db[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
seven_k_same_n_valid$db[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

seven_k_same_n_valid$ascm[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
seven_k_same_n_valid$ascm[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi



### A3 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,100 )
  sigma <- matrix(c(100,0,0,100), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( -150,50 )
  sigma <- matrix(c(80,0,0,80), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( -50,0 )
  sigma <- matrix(c(60,0,0,60), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 200,0 )
  sigma <- matrix(c(40,0,0,40), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(70,0,0,70), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # cluster 6
  n <- 50
  center <- c( 200,100 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data6 <- generateGaussianData(n, center, sigma, 6)
  # cluster 7
  n <- 50
  center <- c( -150,-50 )
  sigma <- matrix(c(90,0,0,90), nrow = 2)
  data7 <- generateGaussianData(n, center, sigma, 7)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5, data6,data7)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)


kmsil_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 10, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 7, nstart = 25, graph = F)



seven_k_same_n_valid$silhouette[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
seven_k_same_n_valid$silhouette[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

seven_k_same_n_valid$ch[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
seven_k_same_n_valid$ch[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

seven_k_same_n_valid$db[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
seven_k_same_n_valid$db[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

seven_k_same_n_valid$ascm[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
seven_k_same_n_valid$ascm[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi


##############
### Case B ###
##############

### B1 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,75 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( -150,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( -50,25 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 200,25 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # cluster 6
  n <- 50
  center <- c( 200,75 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data6 <- generateGaussianData(n, center, sigma, 6)
  # cluster 7
  n <- 50
  center <- c( -150,0 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data7 <- generateGaussianData(n, center, sigma, 7)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5, data6,data7)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)



kmsil_data <- eclust(df, "kmeans", k = 9, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 10, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 9, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 7, nstart = 25, graph = F)



seven_k_same_n_valid$silhouette[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
seven_k_same_n_valid$silhouette[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

seven_k_same_n_valid$ch[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
seven_k_same_n_valid$ch[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

seven_k_same_n_valid$db[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
seven_k_same_n_valid$db[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

seven_k_same_n_valid$ascm[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
seven_k_same_n_valid$ascm[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi



### B2 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,75 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( -150,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( -50,25 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 200,25 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # cluster 6
  n <- 50
  center <- c( 200,75 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data6 <- generateGaussianData(n, center, sigma, 6)
  # cluster 7
  n <- 50
  center <- c( -150,0 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data7 <- generateGaussianData(n, center, sigma, 7)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5, data6,data7)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)



kmsil_data <- eclust(df, "kmeans", k = 9, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 10, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 9, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 7, nstart = 25, graph = F)



seven_k_same_n_valid$silhouette[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
seven_k_same_n_valid$silhouette[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

seven_k_same_n_valid$ch[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
seven_k_same_n_valid$ch[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

seven_k_same_n_valid$db[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
seven_k_same_n_valid$db[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

seven_k_same_n_valid$ascm[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
seven_k_same_n_valid$ascm[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi



### B3 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,75 )
  sigma <- matrix(c(100,0,0,100), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( -150,50 )
  sigma <- matrix(c(80,0,0,80), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( -50,25 )
  sigma <- matrix(c(60,0,0,60), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 200,25 )
  sigma <- matrix(c(40,0,0,40), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(70,0,0,70), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # cluster 6
  n <- 50
  center <- c( 200,75 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data6 <- generateGaussianData(n, center, sigma, 6)
  # cluster 7
  n <- 50
  center <- c( -150,0 )
  sigma <- matrix(c(90,0,0,90), nrow = 2)
  data7 <- generateGaussianData(n, center, sigma, 7)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5, data6,data7)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)


kmsil_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 10, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 9, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 7, nstart = 25, graph = F)



seven_k_same_n_valid$silhouette[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
seven_k_same_n_valid$silhouette[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

seven_k_same_n_valid$ch[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
seven_k_same_n_valid$ch[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

seven_k_same_n_valid$db[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
seven_k_same_n_valid$db[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

seven_k_same_n_valid$ascm[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
seven_k_same_n_valid$ascm[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi



################################################################################
# //////////////////////////////////////////////////////////////////////////// #
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ###
#### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ####
################################################################################

#############################################
### Validations for 8 clustered data sets ###
#############################################


##############
### Case A ###
##############


### A1 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,100 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( -150,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( -50,0 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 200,0 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # cluster 6
  n <- 50
  center <- c( 200,100 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data6 <- generateGaussianData(n, center, sigma, 6)
  # cluster 7
  n <- 50
  center <- c( -150,-50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data7 <- generateGaussianData(n, center, sigma, 7)
  # cluster 8
  n <- 50
  center <- c( 100,-50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data8 <- generateGaussianData(n, center, sigma, 8)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5, data6, data7, data8)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)



eight_k_same_n_valid <- list()
eight_k_same_n_valid$silhouette <- data.frame()
eight_k_same_n_valid$ch <- data.frame()
eight_k_same_n_valid$db <- data.frame()
eight_k_same_n_valid$ascm <- data.frame()


kmsil_data <- eclust(df, "kmeans", k = 7, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 8, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 8, nstart = 25, graph = F)



eight_k_same_n_valid$silhouette[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
eight_k_same_n_valid$silhouette[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

eight_k_same_n_valid$ch[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
eight_k_same_n_valid$ch[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

eight_k_same_n_valid$db[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
eight_k_same_n_valid$db[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

eight_k_same_n_valid$ascm[1,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
eight_k_same_n_valid$ascm[1,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi


names(eight_k_same_n_valid$silhouette) <- c("rand", "mvi")
names(eight_k_same_n_valid$ch) <- c("rand", "mvi")
names(eight_k_same_n_valid$db) <- c("rand", "mvi")
names(eight_k_same_n_valid$ascm) <- c("rand", "mvi")




### A2 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,100 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( -150,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( -50,0 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 200,0 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # cluster 6
  n <- 50
  center <- c( 200,100 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data6 <- generateGaussianData(n, center, sigma, 6)
  # cluster 7
  n <- 50
  center <- c( -150,-50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data7 <- generateGaussianData(n, center, sigma, 7)
  # cluster 7
  n <- 50
  center <- c( 100,-50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data8 <- generateGaussianData(n, center, sigma, 8)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5, data6, data7, data8)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)


kmsil_data <- eclust(df, "kmeans", k = 7, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 10, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 6, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 8, nstart = 25, graph = F)



eight_k_same_n_valid$silhouette[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
eight_k_same_n_valid$silhouette[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

eight_k_same_n_valid$ch[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
eight_k_same_n_valid$ch[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

eight_k_same_n_valid$db[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
eight_k_same_n_valid$db[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

eight_k_same_n_valid$ascm[2,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
eight_k_same_n_valid$ascm[2,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi


### A3 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,100 )
  sigma <- matrix(c(100,0,0,100), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( -150,50 )
  sigma <- matrix(c(80,0,0,80), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( -50,0 )
  sigma <- matrix(c(60,0,0,60), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 200,0 )
  sigma <- matrix(c(40,0,0,40), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(70,0,0,70), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # cluster 6
  n <- 50
  center <- c( 200,100 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data6 <- generateGaussianData(n, center, sigma, 6)
  # cluster 7
  n <- 50
  center <- c( -150,-50 )
  sigma <- matrix(c(90,0,0,90), nrow = 2)
  data7 <- generateGaussianData(n, center, sigma, 7)
  # cluster 7
  n <- 50
  center <- c( 100,-50 )
  sigma <- matrix(c(110,0,0,110), nrow = 2)
  data8 <- generateGaussianData(n, center, sigma, 8)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5, data6, data7, data8)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)



kmsil_data <- eclust(df, "kmeans", k = 7, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 10, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 8, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 8, nstart = 25, graph = F)



eight_k_same_n_valid$silhouette[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
eight_k_same_n_valid$silhouette[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

eight_k_same_n_valid$ch[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
eight_k_same_n_valid$ch[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

eight_k_same_n_valid$db[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
eight_k_same_n_valid$db[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

eight_k_same_n_valid$ascm[3,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
eight_k_same_n_valid$ascm[3,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi



##############
### Case B ###
##############

### B1 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,75 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( -150,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( -50,25 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 200,25 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # cluster 6
  n <- 50
  center <- c( 200,75 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data6 <- generateGaussianData(n, center, sigma, 6)
  # cluster 7
  n <- 50
  center <- c( -150,0 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data7 <- generateGaussianData(n, center, sigma, 7)
  # cluster 8
  n <- 50
  center <- c( 100,0 )
  sigma <- matrix(c(25,0,0,25), nrow = 2)
  data8 <- generateGaussianData(n, center, sigma, 8)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5, data6,data7, data8)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)



kmsil_data <- eclust(df, "kmeans", k = 7, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 10, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 9, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 8, nstart = 25, graph = F)



eight_k_same_n_valid$silhouette[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
eight_k_same_n_valid$silhouette[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

eight_k_same_n_valid$ch[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
eight_k_same_n_valid$ch[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

eight_k_same_n_valid$db[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
eight_k_same_n_valid$db[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

eight_k_same_n_valid$ascm[4,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
eight_k_same_n_valid$ascm[4,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi


### B2 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,75 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( -150,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( -50,25 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 200,25 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # cluster 6
  n <- 50
  center <- c( 200,75 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data6 <- generateGaussianData(n, center, sigma, 6)
  # cluster 7
  n <- 50
  center <- c( -150,0 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data7 <- generateGaussianData(n, center, sigma, 7)
  # cluster 8
  n <- 50
  center <- c( 100,0 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data8 <- generateGaussianData(n, center, sigma, 8)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5, data6,data7, data8)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)




kmsil_data <- eclust(df, "kmeans", k = 7, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 8, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 9, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 8, nstart = 25, graph = F)



eight_k_same_n_valid$silhouette[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
eight_k_same_n_valid$silhouette[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

eight_k_same_n_valid$ch[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
eight_k_same_n_valid$ch[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

eight_k_same_n_valid$db[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
eight_k_same_n_valid$db[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

eight_k_same_n_valid$ascm[5,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
eight_k_same_n_valid$ascm[5,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi



### B3 ###

dataset1 <- {
  # cluster 1
  n <- 50
  center <- c( -50,75 )
  sigma <- matrix(c(100,0,0,100), nrow = 2)
  data1 <- generateGaussianData(n, center, sigma, 1)
  # cluster 2
  n <- 50
  center <- c( -150,50 )
  sigma <- matrix(c(80,0,0,80), nrow = 2)
  data2 <- generateGaussianData(n, center, sigma, 2)
  # cluster 3
  n <- 50
  center <- c( -50,25 )
  sigma <- matrix(c(60,0,0,60), nrow = 2)
  data3 <- generateGaussianData(n, center, sigma, 3)
  # cluster 4
  n <- 50
  center <- c( 200,25 )
  sigma <- matrix(c(40,0,0,40), nrow = 2)
  data4 <- generateGaussianData(n, center, sigma, 4)
  # cluster 5
  n <- 50
  center <- c( 100,50 )
  sigma <- matrix(c(70,0,0,70), nrow = 2)
  data5 <- generateGaussianData(n, center, sigma, 5)
  # cluster 6
  n <- 50
  center <- c( 200,75 )
  sigma <- matrix(c(50,0,0,50), nrow = 2)
  data6 <- generateGaussianData(n, center, sigma, 6)
  # cluster 7
  n <- 50
  center <- c( -150,0 )
  sigma <- matrix(c(90,0,0,90), nrow = 2)
  data7 <- generateGaussianData(n, center, sigma, 7)
  # cluster 8
  n <- 50
  center <- c( 100,0 )
  sigma <- matrix(c(110,0,0,110), nrow = 2)
  data8 <- generateGaussianData(n, center, sigma, 8)
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5, data6,data7, data8)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)


kmsil_data <- eclust(df, "kmeans", k = 7, nstart = 25, graph = F)
kmch_data <- eclust(df, "kmeans", k = 7, nstart = 25, graph = F)
kmdb_data <- eclust(df, "kmeans", k = 8, nstart = 25, graph = F)
kmascm_data <- eclust(df, "kmeans", k = 8, nstart = 25, graph = F)



eight_k_same_n_valid$silhouette[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmsil_data$cluster))
eight_k_same_n_valid$silhouette[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmsil_data$cluster)$vi

eight_k_same_n_valid$ch[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmch_data$cluster))
eight_k_same_n_valid$ch[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmch_data$cluster)$vi

eight_k_same_n_valid$db[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmdb_data$cluster))
eight_k_same_n_valid$db[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmdb_data$cluster)$vi

eight_k_same_n_valid$ascm[6,1] <- rand.index(as.numeric(dataset1$class), as.numeric(kmascm_data$cluster))
eight_k_same_n_valid$ascm[6,2] <- cluster.stats(d = dist(df, method = "euclidean") ,as.numeric(dataset1$class), kmascm_data$cluster)$vi





################################################################################
################################################################################
################################################################################



all_si <- rbind(two_k_same_n_valid$silhouette, three_k_same_n_valid$silhouette, four_k_same_n_valid$silhouette,
                five_k_same_n_valid$silhouette, six_k_same_n_valid$silhouette, seven_k_same_n_valid$silhouette, eight_k_same_n_valid$silhouette)


all_ch <- rbind(two_k_same_n_valid$ch, three_k_same_n_valid$ch, four_k_same_n_valid$ch,
                five_k_same_n_valid$ch, six_k_same_n_valid$ch, seven_k_same_n_valid$ch, eight_k_same_n_valid$ch)


all_db <- rbind(two_k_same_n_valid$db, three_k_same_n_valid$db, four_k_same_n_valid$db,
                five_k_same_n_valid$db, six_k_same_n_valid$db, seven_k_same_n_valid$db, eight_k_same_n_valid$db)


all_ascm <- rbind(two_k_same_n_valid$ascm, three_k_same_n_valid$ascm, four_k_same_n_valid$ascm,
                five_k_same_n_valid$ascm, six_k_same_n_valid$ascm, seven_k_same_n_valid$ascm, eight_k_same_n_valid$ascm)



averages <- data.frame()

# rand indexes
averages[1,1] <- mean(all_si$rand)
averages[2,1] <- mean(all_ch$rand)
averages[3,1] <- mean(all_db$rand)
averages[4,1] <- mean(all_ascm$rand)


# mvi scores
averages[1,2] <- mean(all_si$mvi)
averages[2,2] <- mean(all_ch$mvi)
averages[3,2] <- mean(all_db$mvi)
averages[4,2] <- mean(all_ascm$mvi)


# success rates

averages[1,3] <- 18/42
averages[2,3] <- 22/42
averages[3,3] <- 16/42
averages[4,3] <- 42/42


names(averages) <- c("rand", "mvi", "success_rate")
row.names(averages) <- c("as", "ch", "db", "ascm")

averages

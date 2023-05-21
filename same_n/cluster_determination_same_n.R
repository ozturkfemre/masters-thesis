library(mvtnorm) # data generation
library(clusterSim) # davies bouldin
library(fpc)# ch
library(magrittr) # pipe
library(dplyr) # mutate
library(factoextra) # silhouette

##########################
### Required Functions ###
##########################

#####################
# Calinski Harabasz #
#####################

# The following function calculates ch score of a data set in range 2 to 10 clusters.

ch <- function(data) {
  ch <- c()
  for (i in 2:10) {
    km <- kmeans(data, i) # perform clustering
    ch[i] <- calinhara(data, # data
                       km$cluster, # cluster assignments
                       cn=max(km$cluster) # total cluster number
    )
  }
  ch <-ch[2:10]
  ch
}

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #

####################
# Davies - Bouldin #
####################

# The following function calculates db score of a data set in range 2 to 10 clusters.

db <- function(data) {
  db <- c()
  for (i in 2:10) {
    km <- kmeans(data, i) # perform clustering
    db[i] <- index.DB(data, 
                      km$cluster, # vector of integers indicating the cluster to which each object is allocated
                      centrotypes="centroids", 
                      p=2, # Euclidean distance
                      q=1 # the average distance of objects in the r-th cluster to the centroid or medoid of the r-th cluster
                      )$DB
  }
  db <- db[2:10]
  db
}

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #


###########################################
### Adjusted Silhouette Based on Median ###
###########################################

#' Advanced Silhouette
#'
#' @param df : a dataframe or tibble
#' @param k : cluster number
#' @param algorithm : clustering algorithm (k-means or k-medoids)
#' @description Calculates Advanced Silhouette coefficient for a given number of cluster
#' @details For each cluster i, Advanced Silhouette Coefficient is defined as follows:
#' a(i): average distance between observations in the cluster i and medoid of the cluster i
#' b(i): average distance between observations in the cluster i and medoid of the nearest cluster of cluster i
#'
#' as(i) = (b(i) - a(i)) / max(b(i), a(i))
#' @return advancedsilhouette() returns a Advanced Silhouette Coefficient
#' @export
#'
#' @examples advancedsilhouette(df, 2, "kmeans")
advancedsilhouette <- function(df, k, algorithm) {
  if (algorithm != "kmeans" && algorithm != "pam") {
    stop("Error: invalid clustering algorithm. Please specify 'kmeans' or 'pam'.")
  }
  # calculate clusters according to the clustering algorithm
  if (algorithm == "kmeans") {
    clusters <- stats::kmeans(df, k, iter.max = 50)$cluster
    df_fun <- df
    df_fun[,(ncol(df)+1)] <- clusters
    observations <- c()
    clustermedoids <- data.frame()
    clusteri <- data.frame()
    for (i in 1:k) {
      clusteri <- df[which(df_fun[,ncol(df_fun)] == i),]
      for (j in 1:ncol(clusteri)) {
        clustermedoids[i,j] <- stats::median(clusteri[,j])
      }
      clusteri <- data.frame()
      observations <- c()
    }
    kcenters <- clustermedoids
  } else if (algorithm == "pam") {
    clusters <- cluster::pam(df, k)$clustering
    kcenters <- cluster::pam(df,k)$medoids
  } else {
    stop("Invalid clustering algorithm. Must be 'kmeans' or 'pam'.")
  }
  
  if (k < 2) {
    stop("cluster number k must be more than 1")
  }
  # calculation of cluster variance
  df_fun <- df
  df_fun[,(ncol(df)+1)] <- clusters
  nobservations <- c()
  observations <- c()
  wgss <- c()
  a <- c()
  # for each cluster, finding the average distance between all of the observations and median of the cluster which they are belong to.
  for (i in 1:k) {
    nobservations[i] <- length(which(i==df_fun[,ncol(df_fun)]))
    observations <- which(i==df_fun[,ncol(df_fun)])
    for (o in observations) {
      wgss[o] <- sqrt(sum((df[o,] - kcenters[i,])^2))
    }
    observations <- c()
    a[i] <- sum(wgss, na.rm = T) / nobservations[i]
    wgss <- c()
  }
  a <- sum(a) / k
  a
  # calculation of cluster separation
  cd <- data.frame()
  #finding the distance between clusters
  for (i in 1:k) {
    for (j in 1:k) {
      if (i == j) {
        cd[i,j] <- Inf
      } else {
        cd[i,j] <- sqrt(sum((kcenters[i,] - kcenters[j,])^2))
      }
    }
  }
  # finding the nearest cluster for each cluster
  nearestcluster <- c()
  for (i in 1:k) {
    nearestcluster[i] <- which.min(cd[i,])
  }
  # for each cluster, finding the average distance between all of the observation in the cluster and the nearest cluster's centroid
  wgss <- c()
  b <- c()
  nobservations <- c()
  observations <- c()
  for (i in 1:k) {
    nobservations[i] <- length(which(i==df_fun[,ncol(df_fun)]))
    observations <- which(i==df_fun[,ncol(df_fun)])
    for (o in observations) {
      wgss[o] <- sqrt(sum((df[o,] - kcenters[nearestcluster[i],])^2))
    }
    b[i] <- sum(wgss, na.rm = T) / nobservations[i]
    observations <- c()
    nobservations <- c()
    wgss <- c()
  }
  b <- sum(b) / k
  b
  # calculation of adjusted silhouette coeffiecient based on median
  d <- append(a,b)
  s <- (b-a) / max(d) # (b-a) / max(d)  # b/a
  s
}

# function to calculate advanced silhouette index in a given range.
#' Advanced Silhouette in a given range of cluster numbers
#'
#' @param df a dfframe
#' @param kmin minimum number of cluster
#' @param kmax maximum number of cluster
#' @param algorithm clustering algortihm (k-means or k-medoids)
#' @description see ?advancedsilhouette for detailed information about Advanced Silhouette.
#'
#' @return k_advancedsilhouette returns a vector of advanced silhouette coefficient in a given range.
#' @export
#'
#' @examples
k_advancedsilhouette <- function(df, kmin = 2, kmax = 10, algorithm) {
  as <- c()
  if (algorithm == "kmeans") {
    for (i in kmin:kmax) {
      as[i] <- advancedsilhouette(df, i, algorithm = "kmeans")
    }
  } else if (algorithm == "pam") {
    for (i in kmin:kmax) {
      as[i] <- advancedsilhouette(df, i, algorithm = "pam")
    }
  }
  as <- as[-1]
  as
}
k_advancedsilhouette(df,algorithm = "kmeans")

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #

##########################
### Average Silhouette ###
##########################


si <- function(data) {
  si <- c()
  for (i in 2:10) {
    km <- eclust(data, "kmeans", k = i, graph = F) # perform clustering
    si[i] <- km$silinfo$avg.width
  }
  si <- si[2:10]
  si
}

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #

#######################
### Data Generation ###
#######################

generateGaussianData <- function(n, center, sigma, label) {
  data = rmvnorm(n, mean = center, sigma = sigma)
  data = data.frame(data)
  data = data %>% mutate(class=factor(label))
  data
}


# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #


#####################################################
### Simulation Study of Data sets with 2 clusters ###
#####################################################

two_k_same_n <- data.frame()

caha <- data.frame()
dabo <- data.frame()
silho <- data.frame()
adjustedsi <- data.frame()


colna <- c("k=2", "k=3", "k=4", "k=5", "k=6", "k=7", "k=8" , "k=9", "k=10")
rowna <- c("case_a1", "case_a2", "case_a3", "case_b1", "case_b2", "case_b3")


##############
### Case A ###
##############

variances <- list(a1 = c( c(25,0,0,25), c(25,0,0,25)), a2 = c( c(50,0,0,50), c(50,0,0,50)), a3 = c( c(50,0,0,50), c(75,0,0,75) ),
                  b1 = c( c(25,0,0,25), c(25,0,0,25)), b2 = c( c(50,0,0,50), c(50,0,0,50)), b3 = c( c(50,0,0,50), c(75,0,0,75) )
                  )

x <- list()
for (i in 1:3) {
  dataset1 <- {
    # cluster 1
    n <- 50
    center <- c( 50, 50 )
    sigma <- matrix(c(variances[[i]][1], variances[[i]][2], variances[[i]][3], variances[[i]][4]), nrow = 2)
    data1 <- generateGaussianData(n, center, sigma, 1)
    # cluster 2
    n <- 50
    center <- c( 0, 0 )
    sigma <- matrix(c(variances[[i]][5], variances[[i]][6], variances[[i]][7], variances[[i]][8]), nrow = 2)
    data2 <- generateGaussianData(n, center, sigma, 2)
    # all data
    data <- bind_rows(data1, data2)
    data$dataset <- "1 - Mixture of Gaussians"
    data
  }
  df <- dataset1[,1:2]
  df <- scale(df)
  df <- as.data.frame(df)
  
  x[[i]] <- df
 
  # cluster number determination
  
  sil <- si(df)
  calin <- ch(df)
  davbol <- db(df)
  ascm <- k_advancedsilhouette(df, algorithm = "kmeans")
  
  silho <- rbind(silho, sil)
  caha <- rbind(caha, calin)
  dabo <- rbind(dabo, davbol)
  adjustedsi <- rbind(adjustedsi, ascm)
  df <- data.frame()
}


write.csv(x[[1]], file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/2k_case_a1.csv" )
write.csv(x[[2]], file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/2k_case_a2.csv" )
write.csv(x[[3]], file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/2k_case_a3.csv" )


##############
### Case B ###
##############

for (i in 1:3) {
  dataset1 <- {
    # cluster 1
    n <- 50
    center <- c( 50, 50 )
    sigma <- matrix(c(variances[[i+3]][1], variances[[i+3]][2], variances[[i+3]][3], variances[[i+3]][4]), nrow = 2)
    data1 <- generateGaussianData(n, center, sigma, 1)
    # cluster 2
    n <- 50
    center <- c( 25, 25 )
    sigma <- matrix(c(variances[[i+3]][5], variances[[i+3]][6], variances[[i+3]][7], variances[[i+3]][8]), nrow = 2)
    data2 <- generateGaussianData(n, center, sigma, 2)
    # all data
    data <- bind_rows(data1, data2)
    data$dataset <- "1 - Mixture of Gaussians"
    data
  }
  df <- dataset1[,1:2]
  df <- scale(df)
  df <- as.data.frame(df)
  
  x[[i+3]] <- df
  # cluster number determination
  
  sil <- si(df)
  calin <- ch(df)
  davbol <- db(df)
  ascm <- k_advancedsilhouette(df, algorithm = "kmeans")
  
  silho <- rbind(silho, sil)
  caha <- rbind(caha, calin)
  dabo <- rbind(dabo, davbol)
  adjustedsi <- rbind(adjustedsi, ascm)
  
}

write.csv(x[[4]], file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/2k_case_b1.csv" )
write.csv(x[[5]], file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/2k_case_b2.csv" )
write.csv(x[[6]], file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/2k_case_b3.csv" )

names(silho) <- colna
row.names(silho) <- rowna

names(caha) <- colna
row.names(caha) <- rowna

names(dabo) <- colna
row.names(dabo) <- rowna

names(adjustedsi) <- colna
row.names(adjustedsi) <- rowna


two_k_same_n <- list(avg_silhouette = silho, calinski_harabasz = caha,
                     davies_bouldin = dabo, ascm = adjustedsi)

two_k_same_n

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #

#####################################################
### Simulation Study of Data sets with 3 clusters ###
#####################################################

three_k_same_n <- data.frame()

caha <- data.frame()
dabo <- data.frame()
silho <- data.frame()
adjustedsi <- data.frame()

##############
### Case A ###
##############

variances <- list(a1 = c( c(25,0,0,25), c(25,0,0,25),  c(25,0,0,25)), a2 = c( c(50,0,0,50), c(50,0,0,50),c(50,0,0,50)), a3 = c( c(100,0,0,100), c(80,0,0,80), c(60,0,0,60) ),
                  b1 = c( c(25,0,0,25), c(25,0,0,25),  c(25,0,0,25)), b2 = c( c(50,0,0,50), c(50,0,0,50),c(50,0,0,50)), b3 = c( c(100,0,0,100), c(80,0,0,80), c(60,0,0,60) ) 
                  )

x <- list()

for (i in 1:3) {
  dataset1 <- {
    # cluster 1
    n <- 50
    center <- c( 50,100 )
    sigma <- matrix(c(variances[[i]][1], variances[[i]][2], variances[[i]][3], variances[[i]][4]), nrow = 2)
    data1 <- generateGaussianData(n, center, sigma, 1)
    # cluster 2
    n <- 50
    center <- c( 25,50 )
    sigma <- matrix(c(variances[[i]][5], variances[[i]][6], variances[[i]][7], variances[[i]][8]), nrow = 2)
    data2 <- generateGaussianData(n, center, sigma, 2)
    # cluster 3
    n <- 50
    center <- c( 50,0 )
    sigma <- matrix(c(variances[[i]][9], variances[[i]][10], variances[[i]][11], variances[[i]][12]), nrow = 2)
    data3 <- generateGaussianData(n, center, sigma, 3)
    # all data
    data <- bind_rows(data1, data2, data3)
    data$dataset <- "1 - Mixture of Gaussians"
    data
  }
  df <- dataset1[,1:2]
  df <- scale(df)
  df <- as.data.frame(df)
  
  x[[i]] <- df
  # cluster number determination
  
  sil <- si(df)
  calin <- ch(df)
  davbol <- db(df)
  ascm <- k_advancedsilhouette(df, algorithm = "kmeans")
  
  silho <- rbind(silho, sil)
  caha <- rbind(caha, calin)
  dabo <- rbind(dabo, davbol)
  adjustedsi <- rbind(adjustedsi, ascm)
  
}


write.csv(x[[1]], file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/3k_case_a1.csv" )
write.csv(x[[2]], file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/3k_case_a2.csv" )
write.csv(x[[3]], file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/3k_case_a3.csv" )

##############
### Case B ###
##############

for (i in 1:3) {
  dataset1 <- {
    # cluster 1
    n <- 50
    center <- c( 75,75 )
    sigma <- matrix(c(variances[[i+3]][1], variances[[i+3]][2], variances[[i+3]][3], variances[[i+3]][4]), nrow = 2)
    data1 <- generateGaussianData(n, center, sigma, 1)
    # cluster 2
    n <- 50
    center <- c( 50,50 )
    sigma <- matrix(c(variances[[i+3]][5], variances[[i+3]][6], variances[[i+3]][7], variances[[i+3]][8]), nrow = 2)
    data2 <- generateGaussianData(n, center, sigma, 2)
    # cluster 3
    n <- 50
    center <- c( 25,75 )
    sigma <- matrix(c(variances[[i+3]][9], variances[[i+3]][10], variances[[i+3]][11], variances[[i+3]][12]), nrow = 2)
    data3 <- generateGaussianData(n, center, sigma, 3)
    # all data
    data <- bind_rows(data1, data2, data3)
    data$dataset <- "1 - Mixture of Gaussians"
    data
  }
  df <- dataset1[,1:2]
  df <- scale(df)
  df <- as.data.frame(df)
  
  x[[i+3]] <- df
  
  # cluster number determination
  
  sil <- si(df)
  calin <- ch(df)
  davbol <- db(df)
  ascm <- k_advancedsilhouette(df, algorithm = "kmeans")
  
  silho <- rbind(silho, sil)
  caha <- rbind(caha, calin)
  dabo <- rbind(dabo, davbol)
  adjustedsi <- rbind(adjustedsi, ascm)
  
}


write.csv(x[[4]], file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/3k_case_b1.csv" )
write.csv(x[[5]], file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/3k_case_b2.csv" )
write.csv(x[[6]], file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/3k_case_b3.csv" )


names(silho) <- colna
row.names(silho) <- rowna

names(caha) <- colna
row.names(caha) <- rowna

names(dabo) <- colna
row.names(dabo) <- rowna

names(adjustedsi) <- colna
row.names(adjustedsi) <- rowna


three_k_same_n <- list(avg_silhouette = silho, calinski_harabasz = caha,
                     davies_bouldin = dabo, ascm = adjustedsi)

three_k_same_n


# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #

#####################################################
### Simulation Study of Data sets with 4 clusters ###
#####################################################



caha <- data.frame()
dabo <- data.frame()
silho <- data.frame()
adjustedsi <- data.frame()

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

write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/4k_case_a1.csv" )
# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)


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
write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/4k_case_a2.csv" )
# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)



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
write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/4k_case_a3.csv" )
# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")


silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)

################################################################################
################################################################################

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

write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/4k_case_b1.csv" )
# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)

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


write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/4k_case_b2.csv" )
# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)

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

write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/4k_case_b3.csv" )
# cluster number determination
  
sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)
  



names(silho) <- colna
row.names(silho) <- rowna

names(caha) <- colna
row.names(caha) <- rowna

names(dabo) <- colna
row.names(dabo) <- rowna

names(adjustedsi) <- colna
row.names(adjustedsi) <- rowna


four_k_same_n <- list(avg_silhouette = silho, calinski_harabasz = caha,
                       davies_bouldin = dabo, ascm = adjustedsi)

four_k_same_n



# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
#                                          *******************************************                              #
#                                              ***********************************                                  #


#####################################################
### Simulation Study of Data sets with 5 clusters ###
#####################################################



caha <- data.frame()
dabo <- data.frame()
silho <- data.frame()
adjustedsi <- data.frame()

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


write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/5k_case_a1.csv" )
# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")
ascm
silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)


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


write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/5k_case_a2.csv" )

# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)





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


write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/5k_case_a3.csv" )

# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")
ascm
silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)


################################################################################
################################################################################

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


write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/5k_case_b1.csv" )

# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")
ascm
silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)


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

write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/5k_case_b2.csv" )


# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")
ascm
silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)


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


write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/5k_case_b3.csv" )


# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")
ascm

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)



names(silho) <- colna
row.names(silho) <- rowna

names(caha) <- colna
row.names(caha) <- rowna

names(dabo) <- colna
row.names(dabo) <- rowna

names(adjustedsi) <- colna
row.names(adjustedsi) <- rowna


five_k_same_n <- list(avg_silhouette = silho, calinski_harabasz = caha,
                      davies_bouldin = dabo, ascm = adjustedsi)

five_k_same_n



# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
#                                          *******************************************                              #
#                                              ***********************************                                  #


#####################################################
### Simulation Study of Data sets with 6 clusters ###
#####################################################



caha <- data.frame()
dabo <- data.frame()
silho <- data.frame()
adjustedsi <- data.frame()

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


write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/6k_case_a1.csv" )
# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")
ascm
silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)


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


write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/6k_case_a2.csv" )

# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)




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


write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/6k_case_a3.csv" )

# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)



################################################################################
################################################################################


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
  # all data
  data <- bind_rows(data1, data2, data3, data4, data5, data6)
  data$dataset <- "1 - Mixture of Gaussians"
  data
}
df <- dataset1[,1:2]
df <- scale(df)
df <- as.data.frame(df)


write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/6k_case_b1.csv" )

# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")
ascm
silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)





### A2 ###

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


write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/6k_case_b2.csv" )

# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)





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


write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/6k_case_b3.csv" )

# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)



names(silho) <- colna
row.names(silho) <- rowna

names(caha) <- colna
row.names(caha) <- rowna

names(dabo) <- colna
row.names(dabo) <- rowna

names(adjustedsi) <- colna
row.names(adjustedsi) <- rowna


six_k_same_n <- list(avg_silhouette = silho, calinski_harabasz = caha,
                      davies_bouldin = dabo, ascm = adjustedsi)

six_k_same_n



# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
#                                          *******************************************                              #
#                                              ***********************************                                  #


#####################################################
### Simulation Study of Data sets with 7 clusters ###
#####################################################



caha <- data.frame()
dabo <- data.frame()
silho <- data.frame()
adjustedsi <- data.frame()

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

write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/7k_case_a1.csv" )


# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")
ascm
silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)


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

write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/7k_case_a2.csv" )


# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)




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

write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/7k_case_a3.csv" )



# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")
ascm
silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)



################################################################################
################################################################################


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


write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/7k_case_b1.csv" )


# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)





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


write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/7k_case_b2.csv" )


# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)





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


write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/7k_case_b3.csv" )


# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)



names(silho) <- colna
row.names(silho) <- rowna

names(caha) <- colna
row.names(caha) <- rowna

names(dabo) <- colna
row.names(dabo) <- rowna

names(adjustedsi) <- colna
row.names(adjustedsi) <- rowna


seven_k_same_n <- list(avg_silhouette = silho, calinski_harabasz = caha,
                     davies_bouldin = dabo, ascm = adjustedsi)

seven_k_same_n







# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
#                                          *******************************************                              #
#                                              ***********************************                                  #


#####################################################
### Simulation Study of Data sets with 8 clusters ###
#####################################################



caha <- data.frame()
dabo <- data.frame()
silho <- data.frame()
adjustedsi <- data.frame()

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
  # cluster 7
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

write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/8k_case_a1.csv" )


# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)


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

write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/8k_case_a2.csv" )


# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)




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

write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/8k_case_a3.csv" )


# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)



################################################################################
################################################################################


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


write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/8k_case_b1.csv" )


# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)





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

write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/8k_case_b2.csv" )


# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)





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

write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/8k_case_b3.csv" )


# cluster number determination

sil <- si(df)
calin <- ch(df)
davbol <- db(df)
ascm <- k_advancedsilhouette(df, algorithm = "kmeans")

silho <- rbind(silho, sil)
caha <- rbind(caha, calin)
dabo <- rbind(dabo, davbol)
adjustedsi <- rbind(adjustedsi, ascm)

write.csv(df, file = "C:/Users/Dell/Desktop/R/thesis/same_n/data_sets/8k_case_b3.csv" )

names(silho) <- colna
row.names(silho) <- rowna

names(caha) <- colna
row.names(caha) <- rowna

names(dabo) <- colna
row.names(dabo) <- rowna

names(adjustedsi) <- colna
row.names(adjustedsi) <- rowna


eight_k_same_n <- list(avg_silhouette = silho, calinski_harabasz = caha,
                       davies_bouldin = dabo, ascm = adjustedsi)

eight_k_same_n





################################################################################

library("xlsx")

getwd()
write.xlsx(two_k_same_n, file = "C:/Users/Dell/Desktop/R/thesis/same_n/methods_coefficients/two_k_same_n.xlsx",
           sheetName = "2k", append = FALSE)
write.xlsx(three_k_same_n, file = "C:/Users/Dell/Desktop/R/thesis/same_n/methods_coefficients/three_k_same_n.xlsx",
           sheetName = "3k", append = FALSE)
write.xlsx(four_k_same_n, file = "C:/Users/Dell/Desktop/R/thesis/same_n/methods_coefficients/four_k_same_n.xlsx",
           sheetName = "4k", append = FALSE)
write.xlsx(five_k_same_n, file = "C:/Users/Dell/Desktop/R/thesis/same_n/methods_coefficients/five_k_same_n.xlsx",
           sheetName = "5k", append = FALSE)
write.xlsx(six_k_same_n, file = "C:/Users/Dell/Desktop/R/thesis/same_n/methods_coefficients/six_k_same_n.xlsx",
           sheetName = "6k", append = FALSE)
write.xlsx(seven_k_same_n, file = "C:/Users/Dell/Desktop/R/thesis/same_n/methods_coefficients/seven_k_same_n.xlsx",
           sheetName = "7k", append = FALSE)
write.xlsx(eight_k_same_n, file = "C:/Users/Dell/Desktop/R/thesis/same_n/methods_coefficients/eight_k_same_n.xlsx",
           sheetName = "8k", append = FALSE)






library(clusterSim) # davies bouldin
library(fpc)# ch
library(magrittr) # pipe
library(dplyr) # mutate
library(factoextra) # silhouette
library(readr)


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

#####################
### wine data set ###
#####################


wine <- read_csv("real_data_sets/data_sets/wine.data", 
                 col_names = FALSE)

corrplot::corrplot(cor(wine[-1]), method = "number") # no need for pca


df1 <- scale(wine[-1])


wine_coefficients <- data.frame()


##########################
### average silhouette ###
##########################

wine_coefficients[1:9,1] <-  si(df1)

#########################
### calinski-harabasz ###
#########################

wine_coefficients[1:9,2] <-  ch(df1)

######################
### davies-bouldin ###
######################

wine_coefficients[1:9,3] <-  db(df1)

############
### ascm ###
############

ascm <- k_advancedsilhouette(as.data.frame(df1), algorithm = "kmeans")


wine_coefficients[1:9,4] <-  ascm


# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #

#################
### column_3c ###
#################

column_3C <- read_table("real_data_sets/data_sets/column_3C.dat", 
                        col_names = FALSE)


corrplot::corrplot(cor(column_3C[-7]), method = "number") # pca is necessary


data.pca <- prcomp(column_3C[-7], center = TRUE, scale. = TRUE)
data.pca # 2 principle components are ok
df2 <- predict(data.pca)[,1:2]


column_coefficients <- data.frame()


##########################
### average silhouette ###
##########################

column_coefficients[1:9,1] <-  si(df2)

#########################
### calinski-harabasz ###
#########################

column_coefficients[1:9,2] <-  ch(df2)

######################
### davies-bouldin ###
######################

column_coefficients[1:9,3] <-  db(df2)

############
### ascm ###
############

ascm <- k_advancedsilhouette(as.data.frame(df2), algorithm = "kmeans")


column_coefficients[1:9,4] <-  ascm




# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #

##############
### e.coli ###
##############

ecoli <- read_table("real_data_sets/data_sets/ecoli.data", 
                    col_names = FALSE)

corrplot::corrplot(cor(ecoli[2:8]), method = "number") # pca is not necessary

df3 <- scale(ecoli[2:8])


ecoli_coefficients <- data.frame()


##########################
### average silhouette ###
##########################

ecoli_coefficients[1:9,1] <-  si(df3)

#########################
### calinski-harabasz ###
#########################

ecoli_coefficients[1:9,2] <-  ch(df3)

######################
### davies-bouldin ###
######################

ecoli_coefficients[1:9,3] <-  db(df3)

############
### ascm ###
############

ascm <- k_advancedsilhouette(as.data.frame(df3), algorithm = "kmeans")

ecoli_coefficients[1:9,4] <-  ascm



# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #

############
### iris ###
############

corrplot::corrplot(cor(iris[1:4]), method = "number") # pca is not necessary


data.pca <- prcomp(iris[1:4], center = TRUE, scale. = TRUE)
summary(data.pca) # 2 principle components are ok
df4 <- predict(data.pca)[,1:2]




iris_coefficients <- data.frame()


##########################
### average silhouette ###
##########################

iris_coefficients[1:9,1] <-  si(df4)

#########################
### calinski-harabasz ###
#########################

iris_coefficients[1:9,2] <-  ch(df4)

######################
### davies-bouldin ###
######################

iris_coefficients[1:9,3] <-  db(df4)

############
### ascm ###
############

ascm <- k_advancedsilhouette(as.data.frame(df4), algorithm = "kmeans")
ascm
iris_coefficients[1:9,4] <-  ascm


# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #

################
### haberman ###
################

haberman <- read_csv("real_data_sets/data_sets/haberman.data", 
                     col_names = FALSE)

corrplot::corrplot(cor(haberman[1:3]), method = "number") # pca is not necessary


df5 <- scale(haberman[1:3])



haberman_coefficients <- data.frame()


##########################
### average silhouette ###
##########################

haberman_coefficients[1:9,1] <-  si(df5)

#########################
### calinski-harabasz ###
#########################

haberman_coefficients[1:9,2] <-  ch(df5)

######################
### davies-bouldin ###
######################

haberman_coefficients[1:9,3] <-  db(df5)

############
### ascm ###
############

ascm <- k_advancedsilhouette(as.data.frame(df5), algorithm = "kmeans")
ascm
haberman_coefficients[1:9,4] <-  ascm


# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #

############
### wdbc ###
############

wdbc <- read_csv("real_data_sets/data_sets/wdbc.data", 
                 col_names = FALSE)



corrplot::corrplot(cor(wdbc[3:12]), method = "number") # pca is necessary

data.pca <- prcomp(wdbc[3:12], center = TRUE, scale. = TRUE) # 2 pc is ok
df6 <- predict(data.pca)[,1:2]


wdbc_coefficients <- data.frame()


##########################
### average silhouette ###
##########################

wdbc_coefficients[1:9,1] <-  si(df6)

#########################
### calinski-harabasz ###
#########################

wdbc_coefficients[1:9,2] <-  ch(df6)

######################
### davies-bouldin ###
######################

wdbc_coefficients[1:9,3] <-  db(df6)

############
### ascm ###
############

ascm <- k_advancedsilhouette(as.data.frame(df6), algorithm = "kmeans")
ascm
wdbc_coefficients[1:9,4] <-  ascm



# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #

#####################
### breast tissue ###
#####################

library(readxl)
breast_tissue <- read_excel("real_data_sets/data_sets/BreastTissue.xls", 
                            sheet = "Data")


corrplot::corrplot(cor(breast_tissue[-c(1,2)]), method = "number") # pca is necessary

data.pca <- prcomp(breast_tissue[-c(1,2)], center = TRUE, scale. = TRUE) # 2 pc is ok

df7 <- predict(data.pca)[,1:2]


breast_coefficients <- data.frame()


##########################
### average silhouette ###
##########################

breast_coefficients[1:9,1] <-  si(df7)

#########################
### calinski-harabasz ###
#########################

breast_coefficients[1:9,2] <-  ch(df7)

######################
### davies-bouldin ###
######################

breast_coefficients[1:9,3] <-  db(df7)

############
### ascm ###
############

ascm <- k_advancedsilhouette(as.data.frame(df7), algorithm = "kmeans")
which.max(ascm)

breast_coefficients[1:9,4] <-  ascm



# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #

#############
### glass ###
#############

glass <- read_csv("real_data_sets/data_sets/glass.data", 
                  col_names = FALSE)


corrplot::corrplot(cor(glass[2:10]), method = "number") # pca is necessary

df8 <- scale(glass[2:10])


glass_coefficients <- data.frame()


##########################
### average silhouette ###
##########################

glass_coefficients[1:9,1] <-  si(df8)

#########################
### calinski-harabasz ###
#########################

glass_coefficients[1:9,2] <-  ch(df8)

######################
### davies-bouldin ###
######################

glass_coefficients[1:9,3] <-  db(df8)

############
### ascm ###
############

ascm <- k_advancedsilhouette(as.data.frame(df8), algorithm = "kmeans")
which.max(ascm)

glass_coefficients[1:9,4] <-  ascm




# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #

#############
### yeast ###
#############


yeast <- read_table("real_data_sets/data_sets/yeast.data", 
                    col_names = FALSE)


corrplot::corrplot(cor(yeast[2:7]), method = "number") # pca is not necessary

df9 <- scale(yeast[2:7])


yeast_coefficients <- data.frame()


##########################
### average silhouette ###
##########################

yeast_coefficients[1:9,1] <-  si(df9)

#########################
### calinski-harabasz ###
#########################

yeast_coefficients[1:9,2] <-  ch(df9)

######################
### davies-bouldin ###
######################

yeast_coefficients[1:9,3] <-  db(df9)

############
### ascm ###
############

ascm <- k_advancedsilhouette(as.data.frame(df9), algorithm = "kmeans")
which.max(ascm)


yeast_coefficients[1:9,4] <-  ascm




#################################################################################
################################################################################
################################################################################


library("xlsx")


write.xlsx(wine_coefficients, file = "C:/Users/Dell/Desktop/R/masters-thesis/real_data_sets/methods_coefficients/wine_coefficients.xlsx",
           sheetName = "1", append = FALSE)

write.xlsx(yeast_coefficients, file = "C:/Users/Dell/Desktop/R/masters-thesis/real_data_sets/methods_coefficients/yeast_coefficients.xlsx",
           sheetName = "1", append = FALSE)

write.xlsx(wdbc_coefficients, file = "C:/Users/Dell/Desktop/R/masters-thesis/real_data_sets/methods_coefficients/wdbc_coefficients.xlsx",
           sheetName = "1", append = FALSE)

write.xlsx(iris_coefficients, file = "C:/Users/Dell/Desktop/R/masters-thesis/real_data_sets/methods_coefficients/iris_coefficients.xlsx",
           sheetName = "1", append = FALSE)

write.xlsx(haberman_coefficients, file = "C:/Users/Dell/Desktop/R/masters-thesis/real_data_sets/methods_coefficients/haberman_coefficients.xlsx",
           sheetName = "1", append = FALSE)

write.xlsx(glass_coefficients, file = "C:/Users/Dell/Desktop/R/masters-thesis/real_data_sets/methods_coefficients/glass_coefficients.xlsx",
           sheetName = "1", append = FALSE)

write.xlsx(ecoli_coefficients, file = "C:/Users/Dell/Desktop/R/masters-thesis/real_data_sets/methods_coefficients/ecoli_coefficients.xlsx",
           sheetName = "1", append = FALSE)

write.xlsx(column_coefficients, file = "C:/Users/Dell/Desktop/R/masters-thesis/real_data_sets/methods_coefficients/column_coefficients.xlsx",
           sheetName = "1", append = FALSE)

write.xlsx(breast_coefficients, file = "C:/Users/Dell/Desktop/R/masters-thesis/real_data_sets/methods_coefficients/breast_coefficients.xlsx",
           sheetName = "1", append = FALSE)

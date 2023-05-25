library(mvtnorm)
library(tidyverse)
library(cluster)
library(NbClust)
library(fpc)
library(factoextra)
library(fossil)
library(readr)
library(readxl)



#####################
### Wine Data Set ###
#####################

wine <- read_csv("real_data_sets/data_sets/wine.data", 
                 col_names = FALSE)


df1 <- scale(wine[-1])



wine_valid <- list()

wine_valid$as <- data.frame()
wine_valid$ch <- data.frame()
wine_valid$db <- data.frame()
wine_valid$ascm <- data.frame()


# Since all methods determined the optimal number as three, they will have the same scores. I cluster tne data only once.


km_data <- eclust(df1, "kmeans", k = 3, nstart = 25, graph = F)


wine_valid$as[1,1] <- rand.index(as.numeric(wine$X1), as.numeric(km_data$cluster))
wine_valid$as[1,2] <- cluster.stats(d = dist(df1, method = "euclidean") ,as.numeric(wine$X1), km_data$cluster)$vi

wine_valid$ch[1,1] <- rand.index(as.numeric(wine$X1), as.numeric(km_data$cluster))
wine_valid$ch[1,2] <- cluster.stats(d = dist(df1, method = "euclidean") ,as.numeric(wine$X1), km_data$cluster)$vi

wine_valid$db[1,1] <- rand.index(as.numeric(wine$X1), as.numeric(km_data$cluster))
wine_valid$db[1,2] <- cluster.stats(d = dist(df1, method = "euclidean") ,as.numeric(wine$X1), km_data$cluster)$vi

wine_valid$ascm[1,1] <- rand.index(as.numeric(wine$X1), as.numeric(km_data$cluster))
wine_valid$ascm[1,2] <- cluster.stats(d = dist(df1, method = "euclidean") ,as.numeric(wine$X1), km_data$cluster)$vi


names(wine_valid$as) <- c("rand", "mvi")
names(wine_valid$ch) <- c("rand", "mvi")
names(wine_valid$db) <- c("rand", "mvi")
names(wine_valid$ascm) <- c("rand", "mvi")


#################
### column_3c ###
#################

column_3C <- read_table("real_data_sets/data_sets/column_3C.dat", 
                        col_names = FALSE)


data.pca <- prcomp(column_3C[-7], center = TRUE, scale. = TRUE)
data.pca # 2 principle components are ok
df2 <- predict(data.pca)[,1:2]

library(magrittr)
column_3C %<>% mutate( X7 = case_when(
  X7 == "DH" ~ 1,
  X7 == "SL" ~ 2,
  X7 == "NO" ~ 3
) )

column_valid <- list()

column_valid$as <- data.frame()
column_valid$ch <- data.frame()
column_valid$db <- data.frame()
column_valid$ascm <- data.frame()


kmas_data <- eclust(df2, "kmeans", k = 2, nstart = 25, graph = F)
kmch_data <- eclust(df2, "kmeans", k = 2, nstart = 25, graph = F)
kmdb_data <- eclust(df2, "kmeans", k = 8, nstart = 25, graph = F)
kmascm_data <- eclust(df2, "kmeans", k = 2, nstart = 25, graph = F)

df2 <- as.data.frame(df2)

column_valid$as[1,1] <- rand.index(as.numeric(column_3C$X7), as.numeric(kmas_data$cluster))
column_valid$as[1,2] <- cluster.stats(d = dist(df2, method = "euclidean") ,as.numeric(column_3C$X7), kmas_data$cluster)$vi

column_valid$ch[1,1] <- rand.index(as.numeric(column_3C$X7), as.numeric(kmch_data$cluster))
column_valid$ch[1,2] <- cluster.stats(d = dist(df2, method = "euclidean") ,as.numeric(column_3C$X7), kmch_data$cluster)$vi

column_valid$db[1,1] <- rand.index(as.numeric(column_3C$X7), as.numeric(kmdb_data$cluster))
column_valid$db[1,2] <- cluster.stats(d = dist(df2, method = "euclidean") ,as.numeric(column_3C$X7), kmdb_data$cluster)$vi

column_valid$ascm[1,1] <- rand.index(as.numeric(column_3C$X7), as.numeric(kmascm_data$cluster))
column_valid$ascm[1,2] <- cluster.stats(d = dist(df2, method = "euclidean") ,as.numeric(column_3C$X7), kmascm_data$cluster)$vi


names(column_valid$as) <- c("rand", "mvi")
names(column_valid$ch) <- c("rand", "mvi")
names(column_valid$db) <- c("rand", "mvi")
names(column_valid$ascm) <- c("rand", "mvi")


##############
### e.coli ###
##############

ecoli <- read_table("real_data_sets/data_sets/ecoli.data", 
                  col_names = FALSE)


df3 <- scale(ecoli[2:8])


colnames(ecoli) <- c("squence_name","mcg","gvh","lip","chg","aac","alm1","alm2","site")

ecoli %<>% mutate( site = case_when(
  site == "cp" ~ 1,
  site == "im" ~ 2,
  site == "imS" ~ 3,
  site == "imL" ~ 4,
  site == "imU" ~ 5,
  site == "om" ~ 6,
  site == "omL" ~ 7,
  site == "pp" ~ 8
) )


ecoli_valid <- list()

ecoli_valid$as <- data.frame()
ecoli_valid$ch <- data.frame()
ecoli_valid$db <- data.frame()
ecoli_valid$ascm <- data.frame()


kmas_data <- eclust(df3, "kmeans", k = 3, nstart = 25, graph = F)
kmch_data <- eclust(df3, "kmeans", k = 9, nstart = 25, graph = F)
kmdb_data <- eclust(df3, "kmeans", k = 8, nstart = 25, graph = F)
kmascm_data <- eclust(df3, "kmeans", k = 8, nstart = 25, graph = F)


ecoli_valid$as[1,1] <- rand.index(as.numeric(ecoli$site), as.numeric(kmas_data$cluster))
ecoli_valid$as[1,2] <- cluster.stats(d = dist(df3, method = "euclidean") ,as.numeric(ecoli$site), kmas_data$cluster)$vi

ecoli_valid$ch[1,1] <- rand.index(as.numeric(ecoli$site), as.numeric(kmch_data$cluster))
ecoli_valid$ch[1,2] <- cluster.stats(d = dist(df3, method = "euclidean") ,as.numeric(ecoli$site), kmch_data$cluster)$vi

ecoli_valid$db[1,1] <- rand.index(as.numeric(ecoli$site), as.numeric(kmdb_data$cluster))
ecoli_valid$db[1,2] <- cluster.stats(d = dist(df3, method = "euclidean") ,as.numeric(ecoli$site), kmdb_data$cluster)$vi

ecoli_valid$ascm[1,1] <- rand.index(as.numeric(ecoli$site), as.numeric(kmascm_data$cluster))
ecoli_valid$ascm[1,2] <- cluster.stats(d = dist(df3, method = "euclidean") ,as.numeric(ecoli$site), kmascm_data$cluster)$vi


names(ecoli_valid$as) <- c("rand", "mvi")
names(ecoli_valid$ch) <- c("rand", "mvi")
names(ecoli_valid$db) <- c("rand", "mvi")
names(ecoli_valid$ascm) <- c("rand", "mvi")



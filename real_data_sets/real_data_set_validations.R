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


############
### iris ###
############

data.pca <- prcomp(iris[1:4], center = TRUE, scale. = TRUE)
summary(data.pca) # 2 principle components are ok
data <- iris

df4 <- predict(data.pca)[,1:2]

data %<>% mutate(Species = case_when(
  Species == "setosa" ~ 1,
  Species == "versicolor" ~ 2,
  Species == "virginica" ~ 3
)) 


iris_valid <- list()

iris_valid$as <- data.frame()
iris_valid$ch <- data.frame()
iris_valid$db <- data.frame()
iris_valid$ascm <- data.frame()


kmas_data <- eclust(df4, "kmeans", k = 2, nstart = 25, graph = F)
kmch_data <- eclust(df4, "kmeans", k =2, nstart = 25, graph = F)
kmdb_data <- eclust(df4, "kmeans", k = 2, nstart = 25, graph = F)
kmascm_data <- eclust(df4, "kmeans", k =2, nstart = 25, graph = F)


iris_valid$as[1,1] <- rand.index(as.numeric(data$Species), 
                                 as.numeric(kmas_data$cluster))
iris_valid$as[1,2] <- cluster.stats(d = dist(df4, method = "euclidean") ,
                                    as.numeric(data$Species), kmas_data$cluster)$vi

iris_valid$ch[1,1] <- rand.index(as.numeric(data$Species), as.numeric(kmch_data$cluster))
iris_valid$ch[1,2] <- cluster.stats(d = dist(df4, method = "euclidean") ,
                                    as.numeric(data$Species), kmch_data$cluster)$vi

iris_valid$db[1,1] <- rand.index(as.numeric(data$Species), 
                                 as.numeric(kmdb_data$cluster))
iris_valid$db[1,2] <- cluster.stats(d = dist(df4, method = "euclidean") ,
                                    as.numeric(data$Species), kmdb_data$cluster)$vi

iris_valid$ascm[1,1] <- rand.index(as.numeric(data$Species), 
                                   as.numeric(kmascm_data$cluster))
iris_valid$ascm[1,2] <- cluster.stats(d = dist(df4, method = "euclidean") ,
                                      as.numeric(data$Species), kmascm_data$cluster)$vi


names(iris_valid$as) <- c("rand", "mvi")
names(iris_valid$ch) <- c("rand", "mvi")
names(iris_valid$db) <- c("rand", "mvi")
names(iris_valid$ascm) <- c("rand", "mvi")


################
### haberman ###
################

haberman <- read_csv("real_data_sets/data_sets/haberman.data", 
                     col_names = FALSE)


df5 <- scale(haberman[1:3])


unique(haberman$X4)

haberman_valid <- list()

haberman_valid$as <- data.frame()
haberman_valid$ch <- data.frame()
haberman_valid$db <- data.frame()
haberman_valid$ascm <- data.frame()


kmas_data <- eclust(df5, "kmeans", k = 5, nstart = 25, graph = F)
kmch_data <- eclust(df5, "kmeans", k =5, nstart = 25, graph = F)
kmdb_data <- eclust(df5, "kmeans", k = 4, nstart = 25, graph = F)
kmascm_data <- eclust(df5, "kmeans", k = 5, nstart = 25, graph = F)


haberman_valid$as[1,1] <- rand.index(as.numeric(haberman$X4), 
                                 as.numeric(kmas_data$cluster))
haberman_valid$as[1,2] <- cluster.stats(d = dist(df5, method = "euclidean") ,
                                    as.numeric(haberman$X4), kmas_data$cluster)$vi

haberman_valid$ch[1,1] <- rand.index(as.numeric(haberman$X4), as.numeric(kmch_data$cluster))
haberman_valid$ch[1,2] <- cluster.stats(d = dist(df5, method = "euclidean") ,
                                    as.numeric(haberman$X4), kmch_data$cluster)$vi

haberman_valid$db[1,1] <- rand.index(as.numeric(haberman$X4), 
                                 as.numeric(kmdb_data$cluster))
haberman_valid$db[1,2] <- cluster.stats(d = dist(df5, method = "euclidean") ,
                                    as.numeric(haberman$X4), kmdb_data$cluster)$vi

haberman_valid$ascm[1,1] <- rand.index(as.numeric(haberman$X4), 
                                   as.numeric(kmascm_data$cluster))
haberman_valid$ascm[1,2] <- cluster.stats(d = dist(df5, method = "euclidean") ,
                                      as.numeric(haberman$X4), kmascm_data$cluster)$vi


names(haberman_valid$as) <- c("rand", "mvi")
names(haberman_valid$ch) <- c("rand", "mvi")
names(haberman_valid$db) <- c("rand", "mvi")
names(haberman_valid$ascm) <- c("rand", "mvi")



############
### wdbc ###
############

wdbc <- read_csv("real_data_sets/data_sets/wdbc.data", 
                 col_names = FALSE)


data.pca <- prcomp(wdbc[3:12], center = TRUE, scale. = TRUE) # 2 pc is ok
df6 <- predict(data.pca)[,1:2]

data <- wdbc[1:12] 
names(data) <- c("ID", "Diagnosis", "radius", "texture", "perimeter", "area", "smoothness", "compactness", "concavity", "concave.points", "symmetry", "fractal.dimension" )

data %<>% mutate( Diagnosis = case_when(
  Diagnosis == "M" ~ 1,
  Diagnosis == "B" ~ 2
) )



wdbc_valid <- list()

wdbc_valid$as <- data.frame()
wdbc_valid$ch <- data.frame()
wdbc_valid$db <- data.frame()
wdbc_valid$ascm <- data.frame()


kmas_data <- eclust(df6, "kmeans", k = 2, nstart = 25, graph = F)
kmch_data <- eclust(df6, "kmeans", k =2, nstart = 25, graph = F)
kmdb_data <- eclust(df6, "kmeans", k = 2, nstart = 25, graph = F)
kmascm_data <- eclust(df6, "kmeans", k = 2, nstart = 25, graph = F)


wdbc_valid$as[1,1] <- rand.index(as.numeric(data$Diagnosis), 
                                     as.numeric(kmas_data$cluster))
wdbc_valid$as[1,2] <- cluster.stats(d = dist(df6, method = "euclidean") ,
                                        as.numeric(data$Diagnosis), kmas_data$cluster)$vi

wdbc_valid$ch[1,1] <- rand.index(as.numeric(data$Diagnosis), as.numeric(kmch_data$cluster))
wdbc_valid$ch[1,2] <- cluster.stats(d = dist(df6, method = "euclidean") ,
                                        as.numeric(data$Diagnosis), kmch_data$cluster)$vi

wdbc_valid$db[1,1] <- rand.index(as.numeric(data$Diagnosis), 
                                     as.numeric(kmdb_data$cluster))
wdbc_valid$db[1,2] <- cluster.stats(d = dist(df6, method = "euclidean") ,
                                        as.numeric(data$Diagnosis), kmdb_data$cluster)$vi

wdbc_valid$ascm[1,1] <- rand.index(as.numeric(data$Diagnosis), 
                                       as.numeric(kmascm_data$cluster))
wdbc_valid$ascm[1,2] <- cluster.stats(d = dist(df6, method = "euclidean") ,
                                          as.numeric(data$Diagnosis), kmascm_data$cluster)$vi


names(wdbc_valid$as) <- c("rand", "mvi")
names(wdbc_valid$ch) <- c("rand", "mvi")
names(wdbc_valid$db) <- c("rand", "mvi")
names(wdbc_valid$ascm) <- c("rand", "mvi")


#####################
### breast tissue ###
#####################


breast_tissue <- read_excel("real_data_sets/data_sets/BreastTissue.xls", 
                            sheet = "Data")
data.pca <- prcomp(breast_tissue[-c(1,2)], center = TRUE, scale. = TRUE) # 2 pc is ok
df7 <- predict(data.pca)[,1:2]

breast_tissue %<>% mutate(Class = case_when(
  Class == "car" ~ 1,
  Class == "fad" ~ 2,
  Class == "mas" ~ 3,
  Class == "gla" ~ 4,
  Class == "con" ~ 5,
  Class == "adi" ~ 6
)) 

breast_valid <- list()

breast_valid$as <- data.frame()
breast_valid$ch <- data.frame()
breast_valid$db <- data.frame()
breast_valid$ascm <- data.frame()


kmas_data <- eclust(df7, "kmeans", k = 2, nstart = 25, graph = F)
kmch_data <- eclust(df7, "kmeans", k =5, nstart = 25, graph = F)
kmdb_data <- eclust(df7, "kmeans", k = 3, nstart = 25, graph = F)
kmascm_data <- eclust(df7, "kmeans", k = 6, nstart = 25, graph = F)


breast_valid$as[1,1] <- rand.index(as.numeric(breast_tissue$Class), 
                                 as.numeric(kmas_data$cluster))
breast_valid$as[1,2] <- cluster.stats(d = dist(df7, method = "euclidean") ,
                                    as.numeric(breast_tissue$Class), kmas_data$cluster)$vi

breast_valid$ch[1,1] <- rand.index(as.numeric(breast_tissue$Class), as.numeric(kmch_data$cluster))
breast_valid$ch[1,2] <- cluster.stats(d = dist(df7, method = "euclidean") ,
                                    as.numeric(breast_tissue$Class), kmch_data$cluster)$vi

breast_valid$db[1,1] <- rand.index(as.numeric(breast_tissue$Class), 
                                 as.numeric(kmdb_data$cluster))
breast_valid$db[1,2] <- cluster.stats(d = dist(df7, method = "euclidean") ,
                                    as.numeric(breast_tissue$Class), kmdb_data$cluster)$vi

breast_valid$ascm[1,1] <- rand.index(as.numeric(breast_tissue$Class), 
                                   as.numeric(kmascm_data$cluster))
breast_valid$ascm[1,2] <- cluster.stats(d = dist(df7, method = "euclidean") ,
                                      as.numeric(breast_tissue$Class), kmascm_data$cluster)$vi


names(breast_valid$as) <- c("rand", "mvi")
names(breast_valid$ch) <- c("rand", "mvi")
names(breast_valid$db) <- c("rand", "mvi")
names(breast_valid$ascm) <- c("rand", "mvi")



#############
### glass ###
#############

glass <- read_csv("real_data_sets/data_sets/glass.data", 
                  col_names = FALSE)

df8 <- scale(glass[2:10])


glass %<>% mutate(X11 = case_when(
  X11 == 1 ~ 1,
  X11 == 2 ~ 2,
  X11 == 3 ~ 3,
  X11 == 5 ~ 4,
  X11 == 6 ~ 5,
  X11 == 7 ~ 6
)) 


glass_valid <- list()

glass_valid$as <- data.frame()
glass_valid$ch <- data.frame()
glass_valid$db <- data.frame()
glass_valid$ascm <- data.frame()


kmas_data <- eclust(df8, "kmeans", k = 2, nstart = 25, graph = F)
kmch_data <- eclust(df8, "kmeans", k =4, nstart = 25, graph = F)
kmdb_data <- eclust(df8, "kmeans", k = 5, nstart = 25, graph = F)
kmascm_data <- eclust(df8, "kmeans", k = 6, nstart = 25, graph = F)


glass_valid$as[1,1] <- rand.index(as.numeric(glass$X11), 
                                   as.numeric(kmas_data$cluster))
glass_valid$as[1,2] <- cluster.stats(d = dist(df8, method = "euclidean") ,
                                      as.numeric(glass$X11), kmas_data$cluster)$vi

glass_valid$ch[1,1] <- rand.index(as.numeric(glass$X11), as.numeric(kmch_data$cluster))
glass_valid$ch[1,2] <- cluster.stats(d = dist(df8, method = "euclidean") ,
                                      as.numeric(glass$X11), kmch_data$cluster)$vi

glass_valid$db[1,1] <- rand.index(as.numeric(glass$X11), 
                                   as.numeric(kmdb_data$cluster))
glass_valid$db[1,2] <- cluster.stats(d = dist(df8, method = "euclidean") ,
                                      as.numeric(glass$X11), kmdb_data$cluster)$vi

glass_valid$ascm[1,1] <- rand.index(as.numeric(glass$X11), 
                                     as.numeric(kmascm_data$cluster))
glass_valid$ascm[1,2] <- cluster.stats(d = dist(df8, method = "euclidean") ,
                                        as.numeric(glass$X11), kmascm_data$cluster)$vi


names(glass_valid$as) <- c("rand", "mvi")
names(glass_valid$ch) <- c("rand", "mvi")
names(glass_valid$db) <- c("rand", "mvi")
names(glass_valid$ascm) <- c("rand", "mvi")



#############
### yeast ###
#############


yeast <- read_table("real_data_sets/data_sets/yeast.data", 
                    col_names = FALSE)


df9 <- scale(yeast[2:7])

unique(yeast$X10)

yeast %<>% mutate( X10 = case_when(
  X10 == "MIT" ~ 1,
  X10 == "NUC" ~ 2,
  X10 == "CYT" ~ 3,
  X10 == "ME1" ~ 4,
  X10 == "EXC" ~ 5,
  X10 == "ME2" ~ 6,
  X10 == "ME3" ~ 7,
  X10 == "VAC" ~ 8,
  X10 == "POX" ~ 9,
  X10 == "ERL" ~ 10
) )


yeast_valid <- list()

yeast_valid$as <- data.frame()
yeast_valid$ch <- data.frame()
yeast_valid$db <- data.frame()
yeast_valid$ascm <- data.frame()


kmas_data <- eclust(df9, "kmeans", k = 4, nstart = 25, graph = F)
kmch_data <- eclust(df9, "kmeans", k =8, nstart = 25, graph = F)
kmdb_data <- eclust(df9, "kmeans", k =7, nstart = 25, graph = F)
kmascm_data <- eclust(df9, "kmeans", k = 8, nstart = 25, graph = F)


yeast_valid$as[1,1] <- rand.index(as.numeric(yeast$X10), 
                                  as.numeric(kmas_data$cluster))
yeast_valid$as[1,2] <- cluster.stats(d = dist(df9, method = "euclidean") ,
                                     as.numeric(yeast$X10), kmas_data$cluster)$vi

yeast_valid$ch[1,1] <- rand.index(as.numeric(yeast$X10), as.numeric(kmch_data$cluster))
yeast_valid$ch[1,2] <- cluster.stats(d = dist(df9, method = "euclidean") ,
                                     as.numeric(yeast$X10), kmch_data$cluster)$vi

yeast_valid$db[1,1] <- rand.index(as.numeric(yeast$X10), 
                                  as.numeric(kmdb_data$cluster))
yeast_valid$db[1,2] <- cluster.stats(d = dist(df9, method = "euclidean") ,
                                     as.numeric(yeast$X10), kmdb_data$cluster)$vi

yeast_valid$ascm[1,1] <- rand.index(as.numeric(yeast$X10), 
                                    as.numeric(kmascm_data$cluster))
yeast_valid$ascm[1,2] <- cluster.stats(d = dist(df9, method = "euclidean") ,
                                       as.numeric(yeast$X10), kmascm_data$cluster)$vi


names(yeast_valid$as) <- c("rand", "mvi")
names(yeast_valid$ch) <- c("rand", "mvi")
names(yeast_valid$db) <- c("rand", "mvi")
names(yeast_valid$ascm) <- c("rand", "mvi")





################################################################################
################################################################################
################################################################################

avg_si <- rbind(yeast_valid$as, wine_valid$as, wdbc_valid$as, haberman_valid$as, glass_valid$as, 
              ecoli_valid$as, column_valid$as, breast_valid$as, iris_valid$as)

calin <- rbind(yeast_valid$ch, wine_valid$ch, wdbc_valid$ch, haberman_valid$ch, glass_valid$ch, 
                ecoli_valid$ch, column_valid$ch, breast_valid$ch, iris_valid$ch)

dabo <- rbind(yeast_valid$db, wine_valid$db, wdbc_valid$db, haberman_valid$db, glass_valid$db, 
               ecoli_valid$db, column_valid$db, breast_valid$db, iris_valid$db)

ascm <- rbind(yeast_valid$ascm, wine_valid$ascm, wdbc_valid$ascm, haberman_valid$ascm, glass_valid$ascm, 
                ecoli_valid$ascm, column_valid$ascm, breast_valid$ascm, iris_valid$ascm)


averages <- data.frame()

# rand indexes
averages[1,1] <- mean(avg_si$rand)
averages[2,1] <- mean(calin$rand)
averages[3,1] <- mean(dabo$rand)
averages[4,1] <- mean(ascm$rand)


# mvi scores
averages[1,2] <- mean(avg_si$mvi)
averages[2,2] <- mean(calin$mvi)
averages[3,2] <- mean(dabo$mvi)
averages[4,2] <- mean(ascm$mvi)


# success rates

averages[1,3] <- 2/9
averages[2,3] <- 2/9
averages[3,3] <- 3/9
averages[4,3] <- 5/9


names(averages) <- c("rand", "mvi", "success_rate")
row.names(averages) <- c("as", "ch", "db", "ascm")

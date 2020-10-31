library(bio3d)
library(ggplot2)
library(cluster)
library(arules)
library(arulesViz)
library(caret)
library(dplyr)
library(factoextra)
library(proxy)
library(tidyr)
library(reshape2)
library(NbClust)
library(FunCluster)
library(caret)
library(naivebayes)


setwd("E:/Documents/IST 707")
set.seed(121804)

#load data
census <- read.csv("data/acs2015_census_tract_data.csv")
UWay <- read.csv("data/United_Way_Data.csv")

#isolate complete case
census1 <- census[complete.cases(census),]

#Subset to only KY and TN and remove unnecessary columns
KY_TN <- subset(census1, census1$State == "Kentucky" | census1$State == "Tennessee")
KY_TN <- KY_TN[,-c(1,15,17)]

#Convert percentage columns to the actually number of people represented
columns <- c(6:11,15:27,30:34) 
for(i in columns)
{
  KY_TN[[i]] <- round(KY_TN$TotalPop * KY_TN[[i]] / 100)
}
#Convert avg measures to population measure
columns_mean <- c(13,14,28)
for (i in columns_mean)
{
  KY_TN[[i]] <- KY_TN$TotalPop * KY_TN[[i]]
}

#Condense by county
KY_TN <- KY_TN %>% group_by(State, County) %>% summarise_all(sum)
KY_TN$Income <- KY_TN$Income/KY_TN$TotalPop
KY_TN$IncomePerCap <- KY_TN$IncomePerCap/KY_TN$TotalPop
KY_TN$MeanCommute <- KY_TN$MeanCommute/KY_TN$TotalPop

#Prepare United Way data for integration
UWay <- UWay[,c(2,7,8)]
colnames(UWay) <- c("Gross","County", "State")
UWay$Gross <- as.numeric(gsub("[$,]","",UWay$Gross))
UWay$State <- ifelse(UWay$State == "KY", "Kentucky", "Tennessee")

#Merge United Way data into Census data by county. Doing this for each individual state then combining because of duplicate County names
KY <- subset(KY_TN, KY_TN$State == "Kentucky")
UWay_KY <- subset(UWay, UWay$State == "Kentucky")
KY <- full_join(KY, UWay_KY, by ="County")
TN <- subset(KY_TN, KY_TN$State =="Tennessee")
UWay_TN <- subset(UWay, UWay$State == "Tennessee")
TN <- full_join(TN, UWay_TN, by ="County")
KY_TN_Uway <- rbind(KY, TN)
names(KY_TN_Uway)[names(KY_TN_Uway) == "State.x"] <- "State"
KY_TN_Uway <- KY_TN_Uway[,-36]
KY_TN_Uway$Gross <- replace_na(KY_TN_Uway$Gross,0)
KY_TN_Uway$UW <- ifelse(KY_TN_Uway$Gross == 0,"No","Yes")
KY_TN_Uway$Gross <- NULL

#Begin conversion to format for ARM
KY_TN_Uway_ARM <- KY_TN_Uway[,-c(1,2)]

#Function to quickly discretize with explanatory labels for specific columns and quantile levels
quantisize <- function(Data, columns, prob)
{
  df <- Data
  for (k in columns)
  {
    steps <- 1/prob
    newCol <- c()
    q <- quantile(Data[[k]], probs=c(seq(prob, prob*(steps-1), prob)))
    Labels <- c()
    for (i in 0:steps)
    {
      if (i== steps) {break}
      if (i==0) {j <- 0} else {j <- t(q)[[i]]}
      if (i==steps-1) {m <- max(Data[[k]])} else {m <- t(q)[[i+1]]}
      Labels[[i+1]] <- paste("Q",i+1,"(",k,") = (",j,",",m,"]", sep="")
    }
    Data[[k]] <- cut(Data[[k]], breaks=c(-1, t(q),max(Data[[k]])), labels= Labels)
    df <- Data
  }
  return(df)
}

# Convert values to quantile ranges
KY_TN_Uway_ARM <- quantisize(KY_TN_Uway_ARM, colnames(KY_TN_Uway_ARM[,-c(9,33)]), 0.25)
KY_TN_Uway_ARM$Pacific <- ifelse(KY_TN_Uway_ARM$Pacific == 0, "Pac = No", "Pac = Yes")

#Write data to csv then read in as transactions
write.csv(KY_TN_Uway_ARM, "E:/Documents/IST 707/Data/KY_TN_Uway_ARM.csv")
ARM <- read.transactions("E:/Documents/IST 707/Data/KY_TN_Uway_ARM.csv", format="basket", sep=",", header=TRUE, rm.duplicates=FALSE)

#Function to determine optimal support and confidence levels
ruleCounts <- function(data, target, min, max)
{
  Scounts <- c()
  for (i in 1:9)
  {
    for (j in 1:9)
    {
      Scounts[[9*(i-1)+j]] <- length(apriori(data, parameter=list(support=10*min-i*min, conf=10*min-j*min,minlen=3,maxlen=5), 
                                             appearance=list(rhs=target, default="lhs"))@lhs)
    }
  }
  counts <- matrix(Scounts, nrow=9, ncol=9)
  counts <- data.frame(counts, row.names=c(seq(max,min,-min)))
  colnames(counts) <- c(seq(max,min,-min))
  return(counts)
}

library(plotly)
# Determine rules that lead to Yes and rules that lead to No
ruleCounts(ARM, "Yes", 0.1, 0.9)
# support of 0.1 and confidence of 0.9 seems the best combo
Yesrules <- apriori(ARM, parameter = list(support=0.1, conf=0.9, minlen=3, maxlen=5), appearance = list(rhs="Yes", default="lhs"))
arules::inspect(head(sort(Yesrules, by="lift"),10))
plot(sort(Yesrules, by="count", decreasing=TRUE)[1:20], method="graph", engine="htmlwidget")

ruleCounts(ARM, "No", 0.1,0.9)
# support of 0.2 and confidence of 0.9 are best here
Norules <- apriori(ARM, parameter = list(support=0.2, conf=0.9, minlen=2,maxlen=4), appearance = list(rhs="No", default="lhs"))
arules::inspect(head(sort(Norules, by="lift"),25))
plot(sort(Norules, by="count", decreasing=TRUE)[1:20], method="graph", engine="htmlwidget")

#Clustering
# Silhouette method
fviz_nbclust(KY_TN_Uway[, -c(1,2,35)], kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# Elbow method
fviz_nbclust(KY_TN_Uway[, -c(1,2,35)], kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Gap statistic
fviz_nbclust(KY_TN_Uway[, -c(1,2,35)], kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")

# 4 Clusters
k4 <- kmeans(KY_TN_Uway[,-c(1,2,35)], centers=4)
KY_TN_Uway$k4Cluster <- k4$cluster #Add cluster IDs to main data
clusplot(KY_TN_Uway[,-c(1,2,35)],k4$cluster)

# create plot 
KY_TN_Melted <- data.frame(melt(KY_TN_Uway, id.vars = c("k4Cluster", "State", "County"), na.rm = TRUE))
KY_TN_Melted$k4Cluster <- as.factor(KY_TN_Melted$k4Cluster)
ggplot(KY_TN_Melted, aes(k4Cluster, value)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free_y") +
  labs(title = "Boxplots of variables by cluster") +
  theme(axis.text.y = element_blank())

# show table
table(KY_TN_Uway$k4Cluster, KY_TN_Uway$UW) # you can see that population is dominant in the difference between the clusters when looking at table and graph simultaneously

dMatrix1 <- dist(KY_TN_Uway[,-c(1,2,35)], method="euclidean")
groups1 <- hclust(dMatrix1, method="ward.D")
plot(groups1, cex=0.9, hang=1)

dMatrix2 <- dist(KY_TN_Uway[,-c(1,2,35)], method="cosine")
groups2 <- hclust(dMatrix2, method="ward.D")
plot(groups2, cex=0.9, hang=1)

dMatrix3 <- dist(KY_TN_Uway[,-c(1,2,35)], method="euclidean")
groups3 <- hclust(dMatrix3, method="complete")
plot(groups3, cex=0.9, hang=1)

dMatrix4 <- dist(KY_TN_Uway[,-c(1,2,35)], method="cosine")
groups4 <- hclust(dMatrix4, method="complete")
plot(groups4, cex=0.9, hang=1)
rect.hclust(groups4,k=4)

nor <- function(x)
{
  (x-min(x))/(max(x)-min(x))
}

KY_TN_Norm <- as.data.frame(lapply(KY_TN_Uway[,3:34], nor))
KY_TN_Norm <- data.frame(KY_TN_Norm, UW=KY_TN_Uway$UW)
plot(hclust(dist(KY_TN_Norm[,-33], method="cosine"), method="complete"), cex=0.9,hang=1)
rect.hclust(groups4,k=4)

fviz_nbclust(KY_TN_Norm[,-33], kmeans, method="wss")

# NAIVE BAYES
NB.df <- KY_TN_Norm[,-33]

# Create Column indicating United Way activity
s_num <- round(nrow(KY_TN_Norm)/5)
samples <- sample(1:nrow(KY_TN_Norm),s_num)


train <- KY_TN_Norm[-samples,]
test <- KY_TN_Norm[samples,]

Test_no_labels <- test[,-33]
Test_labels <- test[,33]
Train_no_labels <- train[,-33]
Train_labels <- train[,33]

# Build NB Model
nb<-naiveBayes(UW~., data = train)

# Predict on test data
nb_pred <- predict(nb, test)

# Print confusion matrix of predictions 
nb_pred_table <- table(Test_labels, nb_pred)
print(knitr::kable(nb_pred_table))

# accuracy of Naive Bayes model
nb_acc <- (sum(diag(nb_pred_table)/sum(nb_pred_table)))
print(paste0("Naive Bayes Accuracy: ", nb_acc))

nb1 <- NaiveBayes(UW~., data = train)
plot(nb1)

# re-do NB with training set comprising 50/50 UW=Yes and UW=No
y_num <- round(nrow(KY_TN_Norm[KY_TN_Norm$UW=="Yes",])/2)
y_samples <- sample(1:nrow(KY_TN_Norm[KY_TN_Norm$UW=="Yes",]),y_num)
n_samples <- sample(1:nrow(KY_TN_Norm[KY_TN_Norm$UW=="No",]),y_num)

y_set <- subset(KY_TN_Norm, KY_TN_Norm$UW == "Yes")
n_set <- subset(KY_TN_Norm, KY_TN_Norm$UW == "No")

y_train <- y_set[y_samples,]
n_train <- n_set[n_samples,]
y_test <- y_set[-y_samples,]
n_test <- n_set[-n_samples,]

yn_train <- rbind(y_train, n_train)
yn_test <- rbind(y_test, n_test)

yn_Test_no_labels <- yn_test[,-33]
yn_Test_labels <- yn_test[,33]
yn_Train_no_labels <- yn_train[,-33]
yn_Train_labels <- yn_train[,33]

yn_nb <- naiveBayes(UW~., data=yn_train)
yn_nb_pred <- predict(yn_nb, yn_test, type="class")
(yn_nb_pred_table <- table(yn_Test_labels, yn_nb_pred))
(sum(diag(yn_nb_pred_table))/sum(yn_nb_pred_table))

# SVM
library(e1071)
SVM_fit_P <- svm(UW ~ ., data=train, kernel="polynomial", cross=0.1, scale=TRUE)
SVM_test_P<- predict(SVM_fit_P, Test_no_labels, type="class")
(PTable <- table(SVM_test_P, Test_labels))
(MR_P<- sum(diag(PTable))/sum(PTable))

yn_SVM_fit_P <- svm(UW ~ ., data=yn_train, kernel="polynomial", cross=0.1, scale=TRUE)
yn_SVM_test_P<- predict(yn_SVM_fit_P, yn_Test_no_labels, type="class")
(yn_PTable <- table(yn_SVM_test_P, yn_Test_labels))
(yn_MR_P<- sum(diag(yn_PTable))/sum(yn_PTable))


SVM_fit_L <- svm(UW ~ ., data=train, kernel="linear",cost=0.1, scale=TRUE)
SVM_test_L <- predict(SVM_fit_L, Test_no_labels, type="class")
(LTable <- table(SVM_test_L, Test_labels))
(MR_L<- sum(diag(LTable))/sum(LTable))

SVM_fit_R <- svm(UW ~ ., data=train, kernel="radial", cost=0.1, scale=TRUE)
SVM_test_R <- predict(SVM_fit_R, Test_no_labels, type="class")
(RTable <- table(SVM_test_R, Test_labels))
(MR_R<- sum(diag(RTable))/sum(RTable))

SVM_fit_S <- svm(UW ~. , data=train, kernel="sigmoid", cost=0.1, scale=TRUE)
SVM_test_S <- predict(SVM_fit_S, Test_no_labels, type="class")
(STable <- table(SVM_test_S, Test_labels)) 
(MR_S<- sum(diag(STable))/sum(STable))

# knn 
library(class)
k <- round(sqrt(nrow(KY_TN_Norm)))
knn_fit <- class::knn(Train_no_labels, Test_no_labels, cl=Train_labels, k=k, prob=TRUE)
(KTable <- table(knn_fit, Test_labels))
(MR_K <- sum(diag(KTable))/sum(KTable))

yn_knn_fit <- class::knn(yn_Train_no_labels, yn_Test_no_labels, cl=yn_Train_labels, k=k, prob=TRUE)
(yn_KTable <- table(yn_knn_fit, yn_Test_labels))
(yn_MR_K <- sum(diag(yn_KTable))/sum(yn_KTable))
                       
                                              
# random forest
library(randomForest)
RF <- randomForest(UW~., data=train)
test_pred_RF <- predict(RF, Test_no_labels)
RFTable <- table(test_pred_RF, Test_labels)
hist(treesize(RF))
varImpPlot(RF)
(MR_RF <- sum(diag(RFTable))/sum(RFTable))

yn_RF <- randomForest(UW~., data=yn_train)
yn_test_pred_RF <- predict(yn_RF, yn_Test_no_labels)
yn_RFTable <- table(yn_test_pred_RF, yn_Test_labels)
hist(treesize(yn_RF))
varImpPlot(yn_RF)
(yn_MR_RF <- sum(diag(yn_RFTable))/sum(yn_RFTable))


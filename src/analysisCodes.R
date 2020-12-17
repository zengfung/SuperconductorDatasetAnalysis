########################################################################
# LIBRARIES
library(readr)
library(infotheo)
library(ComplexHeatmap)
library(circlize)
library(MASS)
library(e1071)
library(ggplot2)
library(gridExtra)
library(randomForest)

########################################################################
# IMPORT DATA + DATA CLEANSING
train <- read_csv("train.csv")
unique <- read_csv("unique_m.csv")
train<- as.data.frame(train)
unique <- as.data.frame(unique)
train$number_of_elements <- as.integer(train$number_of_elements)

# rename variables
shortened_name <- c("NoE", "M.AM", "WM.AM", "G.AM", "WG.AM", "E.AM",
                    "WE.AM", "R.AM", "WR.AM", "S.AM", "WS.AM", 
                    "M.F", "WM.F", "G.F", "WG.F", "E.F",
                    "WE.F", "R.F", "WR.F", "S.F", "WS.F",
                    "M.AR", "WM.AR", "G.AR", "WG.AR", "E.AR",
                    "WE.AR", "R.AR", "WR.AR", "S.AR", "WS.AR",
                    "M.D", "WM.D", "G.D", "WG.D", "E.D",
                    "WE.D", "R.D", "WR.D", "S.D", "WS.D",
                    "M.EA", "WM.EA", "G.EA", "WG.EA", "E.EA",
                    "WE.EA", "R.EA", "WR.EA", "S.EA", "WS.EA",
                    "M.FH", "WM.FH", "G.FH", "WG.FH", "E.FH",
                    "WE.FH", "R.FH", "WR.FH", "S.FH", "WS.FH",
                    "M.TC", "WM.TC", "G.TC", "WG.TC", "E.TC",
                    "WE.TC", "R.TC", "WR.TC", "S.TC", "WS.TC",
                    "M.V", "WM.V", "G.V", "WG.V", "E.V",
                    "WE.V", "R.V", "WR.V", "S.V", "WS.V", "CT")

# obtain optimal bin amount for each variable for later use
nbins <- rep(0,ncol(train))
nbins[1] <- 9

for (col in 2:82) {
  h <- hist(train[,col], breaks = "Sturges", plot = FALSE)
  nbins[col] <- length(h$breaks)-1
}
nbins[78] <- 7

# discretize all data (except number_of_elements)
data_disc <- train
data_disc2 <- train
data_disc3 <- unique[,1:87]
names(data_disc) <- shortened_name
names(data_disc2) <- shortened_name
for (col in 2:82){
  data_disc[,col] <- discretize(train[,col], "equalwidth", nbins[col])
  data_disc2[,col] <- discretize(train[,col], "equalfreq", 9)
}
data_disc2[,82] <- ifelse(unique$critical_temp<50, 1,2)
data_disc3[,87] <- ifelse(unique$critical_temp<50, 1,2)

######################################################################
# HISTOGRAMS
histograms <- vector(mode = "list", length = 0)

histograms$NoE <- ggplot(data = NULL, aes(train[,1], fill = as.factor(data_disc2$CT))) + 
  geom_histogram(bins = nbins[1], col = "black") +
  labs(x = "No. of Elements", y = "") + theme(legend.position = "none")
histograms$CT <- ggplot(data = NULL, aes(train[,82], fill = as.factor(data_disc2$CT))) + 
  geom_histogram(bins = nbins[82], col = "black") +
  labs(x = "Critical Temperature", y = "") + theme(legend.position = "none")

histograms$M.AM <- ggplot(data = NULL, aes(train[,2], fill = as.factor(data_disc2$CT))) + 
  geom_histogram(bins = 10, col = "black") +
  labs(x = "Approx. Normal", y = "") + theme(legend.position = "none",
                                             axis.text.x=element_blank(),
                                             axis.ticks.x=element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank())
histograms$WM.AM <- ggplot(data = NULL, aes(train[,3], fill = as.factor(data_disc2$CT))) + 
  geom_histogram(bins = 10, col = "black") +
  labs(x = "Skewed to Left", y = "") + theme(legend.position = "none",
                                             axis.text.x=element_blank(),
                                             axis.ticks.x=element_blank(),
                                             axis.text.y = element_blank(),
                                             axis.ticks.y = element_blank())
histograms$WE.AM <- ggplot(data = NULL, aes(train[,7], fill = as.factor(data_disc2$CT))) + 
  geom_histogram(bins = 10, col = "black") +
  labs(x = "Skewed to Right", y = "") + theme(legend.position = "none", 
                                        axis.text.x=element_blank(),
                                        axis.ticks.x=element_blank(),
                                        axis.text.y= element_blank(),
                                        axis.ticks.y = element_blank())
histograms$R.F <- ggplot(data = NULL, aes(train[,18], fill = as.factor(data_disc2$CT))) + 
  geom_histogram(bins = 10, col = "black") +
  labs(x = "Messy", y = "") + theme(legend.position = "none",
                                    axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.text.y = element_blank())
histograms$R.V <- ggplot(data = NULL, aes(train[,78], fill = as.factor(data_disc2$CT))) + 
  geom_histogram(bins = 10, col = "black") +
  labs(x = "Messy", y = "") + theme(legend.position = "none",axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank())
histograms$S.V <- ggplot(data = NULL, aes(train[,80], fill = as.factor(data_disc2$CT))) + 
  geom_histogram(bins = 10, col = "black") +
  labs(x = "Std", y = "") + theme(legend.position = "none",
                                  axis.text.x=element_blank(),
                                  axis.ticks.x=element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank())

histograms$CT
histograms$NoE #distribution similar to CT
histograms$M.AM #distribution close to normal
histograms$WE.AM #distribution close to normal but skewed to right
histograms$WM.AM #distribution close to normal but skewed to left
histograms$R.F #'messy' distribution
histograms$R.V #'messy' distribution


#######################################################################
# MUTUAL CONDITIONAL ENTROPY
MCE <- matrix(0,nrow=82,ncol=82,
              dimnames=list(names(data_disc), names(data_disc)))

# mutual conditional entropy heatmap
for (i in 1:81) {
  for (j in (i+1):82) {
    condent1 <- (mutinformation(data_disc[,i], data_disc[,i]) - mutinformation(data_disc[,i],data_disc[,j]))/entropy(data_disc[,i])
    condent2 <- (mutinformation(data_disc[,j],data_disc[,j]) - mutinformation(data_disc[,j],data_disc[,i]))/entropy(data_disc[,j])
    mutcondent <- (condent1+condent2)/2
    MCE[i,j] <- mutcondent
    MCE[j,i] <- mutcondent
  }
}

#######################################################################
# TYPE 1 VARIABLES ANALYSIS
type <- c("M","WM","G","WG","E","WE", "R","WR","S","WS")

ht1 <- Heatmap(MCE[2:11,2:11], name= "train",
               cluster_rows = FALSE, cluster_columns = FALSE,
               row_labels = type, column_labels = type,
               column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9),
               column_title = "Atomic Mass")

ht2 <- Heatmap(MCE[12:21,12:21], name= "train",
               cluster_rows = FALSE, cluster_columns = FALSE,
               row_labels = type, column_labels = type,
               column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9),
               column_title = "Fie")

ht3 <- Heatmap(MCE[22:31,22:31], name= "train",
               cluster_rows = FALSE, cluster_columns = FALSE,
               row_labels = type, column_labels = type,
               column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9),
               column_title = "Atomic Radius")

ht4 <- Heatmap(MCE[32:41,32:41], name= "train",
               cluster_rows = FALSE, cluster_columns = FALSE,
               row_labels = type, column_labels = type,
               column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9),
               column_title = "Density")

ht5 <- Heatmap(MCE[42:51,42:51], name= "train",
               cluster_rows = FALSE, cluster_columns = FALSE,
               row_labels = type, column_labels = type,
               column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9),
               column_title = "Electron Affinity")

ht6 <- Heatmap(MCE[52:61,52:61], name= "train",
               cluster_rows = FALSE, cluster_columns = FALSE,
               row_labels = type, column_labels = type,
               column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9),
               column_title = "Fusion Heat")

ht7 <- Heatmap(MCE[62:71,62:71], name= "train",
               cluster_rows = FALSE, cluster_columns = FALSE,
               row_labels = type, column_labels = type,
               column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9),
               column_title = "Thermal Conductivity")

ht8 <- Heatmap(MCE[72:81,72:81], name= "train",
               cluster_rows = FALSE, cluster_columns = FALSE,
               row_labels = type, column_labels = type,
               column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9),
               column_title = "Valence")

htlist1 <- ht1 + ht2 + ht3 + ht4
htlist2 <- ht5 + ht6 + ht7 + ht8
htlist1
htlist2

#######################################################################
# TYPE 2 VARIABLES ANALYSIS
type2 <- c("NoE","AM","F","AR", "D","EA","FH","TC", "V","CT")

ht11 <- Heatmap(MCE[c(1,2,12,22,32,42,52,62,72,82),
                    c(1,2,12,22,32,42,52,62,72,82)],
                name= "train", cluster_rows = FALSE, cluster_columns = FALSE,
                row_labels = type2, column_labels = type2,
                column_title = "Mean")

ht12 <- Heatmap(MCE[c(1,3,13,23,33,43,53,63,73,82),
                    c(1,3,13,23,33,43,53,63,73,82)],
                name= "train", cluster_rows = FALSE, cluster_columns = FALSE,
                row_labels = type2, column_labels = type2,
                column_title = "Wtd Mean")

ht13 <- Heatmap(MCE[c(1,4,14,24,34,44,54,64,74,82),
                    c(1,4,14,24,34,44,54,64,74,82)],
                name= "train", cluster_rows = FALSE, cluster_columns = FALSE,
                row_labels = type2, column_labels = type2,
                column_title = "Gmean")

ht14 <- Heatmap(MCE[c(1,5,15,25,35,45,55,65,75,82),
                    c(1,5,15,25,35,45,55,65,75,82)],
                name= "train", cluster_rows = FALSE, cluster_columns = FALSE,
                row_labels = type2, column_labels = type2,
                column_title = "Wtd Gmean")

ht15 <- Heatmap(MCE[c(1,6,16,26,36,46,56,66,76,82),
                    c(1,6,16,26,36,46,56,66,76,82)],
                name= "train", cluster_rows = FALSE, cluster_columns = FALSE,
                row_labels = type2, column_labels = type2,
                column_title = "Entropy")

ht16 <- Heatmap(MCE[c(1,7,17,27,37,47,57,67,77,82),
                    c(1,7,17,27,37,47,57,67,77,82)],
                name= "train", cluster_rows = FALSE, cluster_columns = FALSE,
                row_labels = type2, column_labels = type2,
                column_title = "Wtd Entropy")

ht17 <- Heatmap(MCE[c(1,8,18,28,38,48,58,68,78,82),
                    c(1,8,18,28,38,48,58,68,78,82)],
                name= "train", cluster_rows = FALSE, cluster_columns = FALSE,
                row_labels = type2, column_labels = type2,
                column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9),
                column_title = "Range")

ht18 <- Heatmap(MCE[c(1,9,19,29,39,49,59,69,79,82),
                    c(1,9,19,29,39,49,59,69,79,82)],
                name= "train", cluster_rows = FALSE, cluster_columns = FALSE,
                row_labels = type2, column_labels = type2,
                column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9),
                column_title = "Wtd Range")

ht19 <- Heatmap(MCE[c(1,10,20,30,40,50,60,70,80,82),
                    c(1,10,20,30,40,50,60,70,80,82)],
                name= "train", cluster_rows = FALSE, cluster_columns = FALSE,
                row_labels = type2, column_labels = type2,
                column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9),
                column_title = "Std")

ht20 <- Heatmap(MCE[c(1,11,21,31,41,51,61,71,81,82),
                    c(1,11,21,31,41,51,61,71,81,82)],
                name= "train", cluster_rows = FALSE, cluster_columns = FALSE,
                row_labels = type2, column_labels = type2,
                column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9),
                column_title = "Wtd Std")

htlist4 <- ht11 + ht12 + ht13
htlist5 <- ht14 + ht15 + ht16
htlist6 <- ht17 + ht18 + ht19 + ht20
htlist4
htlist5
htlist6

#######################################################################
# CRITICAL TEMPERATURE CLASSIFICATION BY INFORMATION FLOW
col_fun = colorRamp2(c(1, 9), c("yellow", "red"))

# create heatmap annotation
ha_left <- HeatmapAnnotation(CT = as.factor(data_disc2$CT),
                        col = list(CT = c("1" = "green", "2" = "blue")),
                        annotation_name_side = "left",
                        annotation_legend_param = list(
                          CT = list(direction = "horizontal")))

ha_right <- HeatmapAnnotation(CT = as.factor(data_disc2$CT),
                         col = list(CT = c("1" = "green", "2" = "blue")),
                         annotation_name_side = "right",
                         annotation_legend_param = list(
                           CT = list(direction = "horizontal")))


# heatmap from largest absolute values of correlation
ht21 <- Heatmap(t(as.matrix(data_disc2[,c(71,28,7,1,18,76,57,36,46)])),
                name="train", column_split = 3, row_names_side = "left",
                top_annotation=ha_left, col = col_fun, cluster_rows = FALSE,
                column_title = "Highest Correlation",
                row_names_gp = gpar(fontsize = 9),
                column_title_side = "bottom",
                heatmap_legend_param = list(direction = "horizontal"))

# heatmap from smallest values of mutual conditional entropy
ht22 <- Heatmap(t(as.matrix(data_disc2[,c(34,68,18,73,57,7,1,28,48)])),
                name="train", column_split = 2, row_names_side = "right",
                top_annotation=ha_right, col = col_fun, cluster_rows = FALSE,
                column_title = "Smallest Mutual Cond. Entropy",
                row_names_gp = gpar(fontsize = 9),
                column_title_side = "bottom",
                heatmap_legend_param = list(direction = "horizontal"))

htlist7 <- ht21 + ht22

draw(htlist7, merge_legend = TRUE, heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")

#######################################################################
# OVERALL HEATMAP FOR TRAIN
ha2 <- rowAnnotation(CT = train$critical_temp,
                     annotation_legend_param = list(
                       CT = list(direction = "horizontal")))
crit_temp_order <- order(train$critical_temp)

ht31 <- Heatmap(as.matrix(data_disc2[,1:81]), name = "train",
        left_annotation = ha2, row_order = crit_temp_order,
        column_names_gp = gpar(fontsize = 6), column_km = 3,
        heatmap_legend_param = list(direction = "horizontal"))
draw(ht31, merge_legend = TRUE, heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")

#######################################################################
# OVERALL HEATMAP FOR UNIQUE
ha4 <- HeatmapAnnotation(CT = unique$critical_temp,
                         annotation_name_side = "right",
                     annotation_legend_param = list(
                       CT = list(direction = "horizontal")))
compressed_unique <- data.frame(O = unique$O, Hg = unique$Hg,
                           Cu = unique$Cu, Tl = unique$Tl,
                           Ba = unique$Ba, Pb = unique$Pb,
                           Sr = unique$Sr, Nb = unique$Nb,
                           Ca = unique$Ca, Y = unique$Y,
                           La = unique$La, Bi = unique$Bi,
                           Fe = unique$Fe, As = unique$As)

ht41 <- Heatmap(t(as.matrix(compressed_unique)), name = "count",
        top_annotation = ha4, column_order = crit_temp_order,
        row_names_gp = gpar(fontsize = 8), clustering_method_columns = "complete",
        heatmap_legend_param = list(direction = "horizontal"))
draw(ht41, merge_legend = TRUE, heatmap_legend_side = "bottom", 
         annotation_legend_side = "bottom")

#######################################################################
# PREVALENCE OF ELEMENTS
overall_elemcomp <- colSums(unique[,1:86] != 0)
high_elemcomp <- colSums(unique[which(unique[,87] >= 50),1:86] != 0)
low_elemcomp <- colSums(unique[which(unique[,87] < 50),1:86] != 0)

prop_overall_elemcomp <- overall_elemcomp/nrow(unique)*100
prop_high_elemcomp <- high_elemcomp/sum(unique[,87] >= 50)*100
prop_low_elemcomp <- low_elemcomp/sum(unique[,87] < 50)*100

element_proportion <- data.frame(element = names(unique)[1:86],
                                 overall = prop_overall_elemcomp,
                                 high = prop_high_elemcomp,
                                 low = prop_low_elemcomp)
top10_overall <- order(prop_overall_elemcomp, decreasing = TRUE)[1:10]
top10_high <- order(prop_high_elemcomp, decreasing = TRUE)[1:10]
top10_low <- order(prop_low_elemcomp, decreasing = TRUE)[1:10]

ggplot(data = NULL, aes(x = element_proportion$element[top10_overall], 
                        y = element_proportion$overall[top10_overall])) +
  geom_bar(stat = "identity", fill = "steelblue", col = "black") + 
  geom_text(aes(label=round(element_proportion$overall[top10_overall])),
            position=position_nudge(y = 10), size = 3) +
  scale_x_discrete(limits=element_proportion$element[rev(top10_overall)]) +
  scale_y_continuous(limits = c(0,100)) + coord_flip() +
  labs(x = "Top 10 Elements", y = "Proportion (%)")

ggplot(data = NULL, aes(x = element_proportion$element[top10_high], 
                        y = element_proportion$high[top10_high])) +
  geom_bar(stat = "identity", fill = "steelblue", col = "black") + 
  geom_text(aes(label=round(element_proportion$high[top10_high])),
            position=position_nudge(y = 12), size = 3) +
  scale_x_discrete(limits=element_proportion$element[rev(top10_high)]) +
  scale_y_continuous(limits = c(0,120), breaks = c(0,25,50,75,100)) + coord_flip() +
  labs(x = "Top 10 Elements", y = "Proportion (%)")

ggplot(data = NULL, aes(x = element_proportion$element[top10_low], 
                        y = element_proportion$low[top10_low])) +
  geom_bar(stat = "identity", fill = "steelblue", col = "black") + 
  geom_text(aes(label=round(element_proportion$low[top10_low])),
            position=position_nudge(y = 10),size = 3) +
  scale_x_discrete(limits=element_proportion$element[rev(top10_low)]) +
  scale_y_continuous(limits = c(0,100)) + coord_flip() +
  labs(x = "Top 10 Elements", y = "Proportion (%)")


#######################################################################
# DATA CLEANSING FOR CLASSIFICATION
zero_index <- c()
for (col in 1:ncol(unique)) {
  if (all(unique[,col]==0)) {
    zero_index <- c(zero_index,col)
  }
}
nnz_unique <- unique[,-zero_index]

fulldata <- cbind(nnz_unique[,1:77], train, data_disc2$CT)
names(fulldata) <- c(names(nnz_unique)[1:77], names(data_disc2), "CT_disc")
fulldata$CT_disc <- as.factor(fulldata$CT_disc)

# split full data into 3 classes based on CT_disc
class1 <- fulldata[which(fulldata$CT_disc==1),]
class2 <- fulldata[which(fulldata$CT_disc==2),]

# randomly select 80% from each class to form training data
train1_index <- sample(1:nrow(class1),floor(0.8*nrow(class1)))
train2_index <- sample(1:nrow(class2),floor(0.8*nrow(class2)))
data_train <- rbind(class1[train1_index,], class2[train2_index,])

# remaining 20% as test data
data_test <- rbind(class1[-train1_index,], class2[-train2_index,])

#######################################################################
# SUPERVISED LEARNING: LINEAR REGRESSION
# full model
fullmodel <- lm(CT ~ ., data = data_train[,1:159])

# make predictions and classify them into levels
lmPred <- predict(fullmodel, data_test[,1:159])
lmPred_class <- rep(0,length(lmPred))
for (i in 1:length(lmPred)) {
  if (lmPred[i] < 50) {
    lmPred_class[i] = 1
  } else {
    lmPred_class[i] = 2
  }
}

# confusion matrix between actual vs predicted
confmat_lm <- table(true = data_test$CT_disc,pred = lmPred_class)
confmat_lm

# error rate
errorrate_lm <- 1 - sum(diag(confmat_lm))/sum(confmat_lm)
errorrate_lm

#######################################################################
# SUPERVISED LEARNING: LINEAR DISCRIMINANT ANALYSIS (LDA)
ldamodel <- lda(CT_disc ~ ., data_train[,c(1:158,160)])
ldaPred <- predict(ldamodel,data_test[,c(1:158,160)])
confmat_lda <- table(true = data_test$CT_disc, pred = ldaPred$class)
confmat_lda

# error rate
errorrate_lda <- 1 - sum(diag(confmat_lda))/sum(confmat_lda)
errorrate_lda

#######################################################################
# SUPERVISED LEARNING: RANDOM FOREST
rf_fit = randomForest(CT_disc ~ ., data=data_train[,c(1:158,160)], ntree=100, importance=TRUE)

rf_pred = predict(rf_fit, data_test[,c(1:158,160)])
confmat_rf = table(true = data_test$CT_disc, pred = rf_pred)
confmat_rf

errorrate_rf = 1 - sum(diag(confmat_rf))/sum(confmat_rf)
errorrate_rf

library(dplyr)
#library(maftools)
library(caTools)
library(randomForest)
library(performanceEstimation)
library(data.table)

#--------------------------------------------------PROSTATE---------------------------------------------------

infoPRCA <- read.delim("C:/Users/Giulia/Desktop/Bachelor/FBK/Data/prad_tcga/data_clinical_patient.txt", header = FALSE)
genePRCA <- read.delim("C:/Users/Giulia/Desktop/Bachelor/FBK/Data/prad_tcga/data_mrna_seq_v2_rsem.txt", header = FALSE)
methPRCA <- read.delim("C:/Users/Giulia/Desktop/Bachelor/FBK/Data/prad_tcga/data_methylation_hm450.txt", header = FALSE)
cnaPRCA <- read.delim("C:/Users/Giulia/Desktop/Bachelor/FBK/Data/prad_tcga/data_linear_cna.txt", header = FALSE)

#pre-processing patients info
infoPRCA <- infoPRCA %>% .[-c(1, 2, 3, 4), -1] %>% select(c(1, 6, 7, 8))
colnames(infoPRCA) <- infoPRCA[1, ]
infoPRCA <- infoPRCA[-1, ]

#pre-processing genes dataset
genePRCA[1,] <- gsub("-01", "", genePRCA[1,])
genePRCA <- genePRCA %>% t(.) %>% .[-2, ] %>% as.data.frame(.)
genePRCA[1, 1] <- "PATIENT_ID"
colnames(genePRCA) <- genePRCA[1, ]
genePRCA <- genePRCA[-1, ]
resultsP <- read.delim("C:/Users/Giulia/Desktop/Bachelor/FBK/Data/prad_tcga/data_mutsig.txt")
bestGenesP <- as.vector(t(subset(resultsP, (rank <= 30), select = c(gene))))
colnames(genePRCA)[which(names(genePRCA) == "VPS35L")] <- "C16orf62"
colnames(genePRCA)[which(names(genePRCA) == "TEP1")] <- "PTEN"
colnames(genePRCA)[which(names(genePRCA) == "FGS1")] <- "MED12"
colnames(genePRCA)[which(names(genePRCA) == "OR4P3P")] <- "OR4P4"
genePRCA <- select(genePRCA, all_of(append("PATIENT_ID", bestGenesP)))

#pre-processing DNA methylation dataset
methPRCA[1,] <- gsub("-01", "", methPRCA[1,])
methPRCA <- methPRCA %>% t(.) %>% .[-2, ] %>% as.data.frame(.)
methPRCA[1, 1] <- "PATIENT_ID"
colnames(methPRCA) <- methPRCA[1, ]
methPRCA <- methPRCA[-1, ]
resultsP <- read.delim("C:/Users/Giulia/Desktop/Bachelor/FBK/Data/prad_tcga/data_mutsig.txt")
bestGenesP <- as.vector(t(subset(resultsP, (rank <= 30), select = c(gene))))
colnames(methPRCA)[which(names(methPRCA) == "OPA1")] <- "MED12"
bestGenesP <- bestGenesP[-c(2, 13, 18, 19, 22, 23, 26)]
methPRCA <- select(methPRCA, all_of(append("PATIENT_ID", bestGenesP)))

#pre-processing CNA dataset
cnaPRCA[1,] <- gsub("-01", "", cnaPRCA[1,])
cnaPRCA <- cnaPRCA %>% t(.) %>% .[-2, ] %>% as.data.frame(.)
cnaPRCA[1, 1] <- "PATIENT_ID"
colnames(cnaPRCA) <- cnaPRCA[1, ]
cnaPRCA <- cnaPRCA[-1, ]
resultsP <- read.delim("C:/Users/Giulia/Desktop/Bachelor/FBK/Data/prad_tcga/data_mutsig.txt")
bestGenesP <- as.vector(t(subset(resultsP, (rank <= 30), select = c(gene))))
bestGenesP <- bestGenesP[-17]
cnaPRCA <- select(cnaPRCA, all_of(append("PATIENT_ID", bestGenesP)))

#dataset creation for random forest
dataPRCA <- merge(genePRCA, infoPRCA, by = "PATIENT_ID")
dataPRCA <- merge(methPRCA, infoPRCA, by = "PATIENT_ID")
dataPRCA <- merge(cnaPRCA, infoPRCA, by = "PATIENT_ID")

#unified dataset creation
dataPRCA <- merge(genePRCA, infoPRCA, by = "PATIENT_ID")
dataPRCA <- merge(dataPRCA, methPRCA, by = "PATIENT_ID")
dataPRCA <- merge(dataPRCA, cnaPRCA, by = "PATIENT_ID")

dataPRCA <- dataPRCA %>% .[, -1] %>% sapply(., as.numeric) %>% as.data.frame(.)
dataPRCA <- replace(dataPRCA, is.na(dataPRCA), 0)
dataPRCA <- dataPRCA[, !sapply(dataPRCA, function(x) median(x) == 0)]
dataPRCA <- dataPRCA[dataPRCA$GLEASON_SCORE %in% c(7, 9), ]
dataP1 <- copy(dataPRCA)
dataP1[, 'GLEASON_SCORE'] <- as.factor(dataP1[, 'GLEASON_SCORE'])
dataP1 <- dataP1[, -c(28, 29)]
dataP2 <- dataPRCA[dataPRCA$GLEASON_SCORE %in% c(7), ]
dataP2[, 'GLEASON_PATTERN_PRIMARY'] <- as.factor(dataP2[, 'GLEASON_PATTERN_PRIMARY'])
dataP2 <- dataP2[, -c(29, 30)]

#smote to balance the number of samples
dataP1 <- smote(GLEASON_SCORE ~., dataP1, perc.over = 7, k = 5, perc.under = 1.2)
dataP2 <- smote(GLEASON_PATTERN_PRIMARY ~., dataP2, perc.over = 8, k = 5, perc.under = 1.25)

#Random Forest
#Gleason Score 7 vs 9
split1 <- sample.split(dataP1, SplitRatio = 0.7)
train1 <- subset(dataP1, split1 == "TRUE")
test1 <- subset(dataP1, split1 == "FALSE")
set.seed(120)  # Setting seed
rfP1 = randomForest(x = train1[-28], y = train1$GLEASON_SCORE, ntree = 500)
y_pred = predict(rfP1, newdata = test1[-28])
confusion_mtx = caret::confusionMatrix(test1[, 28], y_pred)
confusion_mtx
plot(rfP1)
#Gleason Pattern 3+4 vs 4+3
split2 <- sample.split(dataP2, SplitRatio = 0.7)
train2 <- subset(dataP2, split2 == "TRUE")
test2 <- subset(dataP2, split2 == "FALSE")
set.seed(120)  # Setting seed
rfP2 = randomForest(x = train2[-28], y = train2$GLEASON_PATTERN_PRIMARY, ntree = 500)
y_pred = predict(rfP2, newdata = test2[-28])
confusion_mtx = caret::confusionMatrix(test2[, 28], y_pred)
confusion_mtx
plot(rfP2)

library(dplyr)
library(maftools)
library(caTools)
library(randomForest)
library(performanceEstimation)
library(data.table)

#-----------------------------------------------------BREAST--------------------------------------------------------

infoBRCA <- read.delim("C:/Users/Giulia/Desktop/Bachelor/FBK/Data/brca_tcga_pub2015/data_clinical_patient.txt", header = FALSE)
geneBRCA <- read.delim("C:/Users/Giulia/Desktop/Bachelor/FBK/Data/brca_tcga_pub2015/data_mrna_seq_v2_rsem.txt", header = FALSE)
methBRCA <- read.delim("C:/Users/Giulia/Desktop/Bachelor/FBK/Data/brca_tcga_pub2015/data_methylation_hm450.txt", header = FALSE)
cnaBRCA <- read.delim("C:/Users/Giulia/Desktop/Bachelor/FBK/Data/brca_tcga_pub2015/data_linear_cna.txt", header = FALSE)

#pre-processing patients info
infoBRCA <- infoBRCA %>% .[-c(1, 2, 3, 4), ] %>% select(c(1, 33))
colnames(infoBRCA) <- infoBRCA[1, ]
infoBRCA <- infoBRCA %>% .[-1, ] %>% 
  rename(TUMOR_STAGE  = AJCC_PATHOLOGIC_TUMOR_STAGE)

#pre-processing genes dataset
geneBRCA[1,] <- gsub("-01", "", geneBRCA[1,])
geneBRCA <- geneBRCA %>% t(.) %>% .[-2, ] %>% as.data.frame(.)
geneBRCA[1, 1] <- "PATIENT_ID"
colnames(geneBRCA) <- geneBRCA[1, ]
geneBRCA <- geneBRCA[-1, ]
resultsB <- read.delim("C:/Users/Giulia/Desktop/Bachelor/FBK/Data/brca_tcga_pub2015/data_mutsig.txt")
bestGenesB <- as.vector(t(subset(resultsB, (rank <= 30), select = c(gene))))
colnames(geneBRCA)[which(names(geneBRCA) == "KMT2C")] <- "MLL3"
colnames(geneBRCA)[which(names(geneBRCA) == "MFRP")] <- "C1QTNF5"
geneBRCA <- select(geneBRCA, all_of(append("PATIENT_ID", bestGenesB)))

#pre-processing DNA methylation dataset
methBRCA[1,] <- gsub("-01", "", methBRCA[1,])
methBRCA <- methBRCA %>% t(.) %>% .[-2, ] %>% as.data.frame(.)
methBRCA[1, 1] <- "PATIENT_ID"
colnames(methBRCA) <- methBRCA[1, ]
methBRCA <- methBRCA[-1, ]
resultsB <- read.delim("C:/Users/Giulia/Desktop/Bachelor/FBK/Data/brca_tcga_pub2015/data_mutsig.txt")
bestGenesB <- as.vector(t(subset(resultsB, (rank <= 30), select = c(gene))))
bestGenesB <- bestGenesB[-c(22, 24)]
colnames(methBRCA)[which(names(methBRCA) == "MFRP")] <- "C1QTNF5"
methBRCA <- select(methBRCA, all_of(append("PATIENT_ID", bestGenesB)))

#pre-processing CNA dataset
cnaBRCA[1,] <- gsub("-01", "", cnaBRCA[1,])
cnaBRCA <- cnaBRCA %>% t(.) %>% .[-2, ] %>% as.data.frame(.)
cnaBRCA[1, 1] <- "PATIENT_ID"
colnames(cnaBRCA) <- cnaBRCA[1, ]
cnaBRCA <- cnaBRCA[-1, ]
resultsB <- read.delim("C:/Users/Giulia/Desktop/Bachelor/FBK/Data/brca_tcga_pub2015/data_mutsig.txt")
bestGenesB <- as.vector(t(subset(resultsB, (rank <= 30), select = c(gene))))
colnames(cnaBRCA)[which(names(cnaBRCA) == "KMT2C")] <- "MLL3"
cnaBRCA <- select(cnaBRCA, all_of(append("PATIENT_ID", bestGenesB)))

#data frame creation for random forest
dataBRCA <- merge(geneBRCA, infoBRCA, by = "PATIENT_ID")
dataBRCA <- merge(methBRCA, infoBRCA, by = "PATIENT_ID")
dataBRCA <- merge(cnaBRCA, infoBRCA, by = "PATIENT_ID")

#unified data frame creation
dataBRCA <- merge(geneBRCA, infoBRCA, by = "PATIENT_ID")
dataBRCA <- merge(dataBRCA, methBRCA, by = "PATIENT_ID")
dataBRCA <- merge(dataBRCA, cnaBRCA, by = "PATIENT_ID")

dataBRCA <- dataBRCA[dataBRCA$TUMOR_STAGE %in% c("Stage IIA", "Stage IIB", "Stage IIIA"), ]
dataBRCA <- dataBRCA %>% 
  mutate(TUMOR_STAGE = gsub("Stage IIA", "1", TUMOR_STAGE)) %>%
  mutate(TUMOR_STAGE = gsub("Stage IIB", "2", TUMOR_STAGE)) %>% 
  mutate(TUMOR_STAGE = gsub("Stage IIIA", "3", TUMOR_STAGE))
dataBRCA <- dataBRCA %>% .[, -1] %>% sapply(., as.numeric) %>% as.data.frame(.)
dataBRCA <- replace(dataBRCA, is.na(dataBRCA), 0)
dataBRCA <- dataBRCA[, !sapply(dataBRCA, function(x) median(x) == 0)]
dataB1 <- copy(dataBRCA)
dataB1[, 29] <- gsub("2", "1", dataB1[, 29])
dataB1[, 'TUMOR_STAGE'] <- as.factor(dataB1[, 'TUMOR_STAGE'])
dataB2 <- dataBRCA[dataBRCA$TUMOR_STAGE %in% c(1, 2), ]
dataB2[, 'TUMOR_STAGE'] <- as.factor(dataB2[, 'TUMOR_STAGE'])

#smote to balance the samples number
dataB1 <- smote(TUMOR_STAGE ~., dataB1, perc.over = 8, k = 5, perc.under = 1.2)
dataB2 <- smote(TUMOR_STAGE ~., dataB2, perc.over = 4, k = 5, perc.under = 1.25)

#Random Forest
#Tumor Stage II vs III
split <- sample.split(dataB1, SplitRatio = 0.7)
train <- subset(dataB1, split == "TRUE")
test <- subset(dataB1, split == "FALSE")
set.seed(120)  # Setting seed
rfB1 = randomForest(x = train[-29], y = train$TUMOR_STAGE, ntree = 500)
y_pred = predict(rfB1, newdata = test[-29])
confusion_mtx = caret::confusionMatrix(test[, 29], y_pred)
confusion_mtx
plot(rfB1)
#Tumor Stage IIA vs IIB
split <- sample.split(dataB2, SplitRatio = 0.7)
train <- subset(dataB2, split == "TRUE")
test <- subset(dataB2, split == "FALSE")
set.seed(120)  # Setting seed
rfB2 = randomForest(x = train[-29], y = train$TUMOR_STAGE, ntree = 500)
y_pred = predict(rfB2, newdata = test[-29])
confusion_mtx = caret::confusionMatrix(test[, 29], y_pred)
confusion_mtx
plot(rfB2)


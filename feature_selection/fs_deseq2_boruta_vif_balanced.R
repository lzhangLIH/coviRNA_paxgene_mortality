
# packages ----------------------------------------------------------------

# BiocManager::install('DESeq2')
library(DESeq2)
# library(edgeR)
# install.packages('Boruta')
library(Boruta)
# install.packages('openxlsx')
library(openxlsx)
# install.packages('car')
library(car)
# install.packages('imbalance')
# library(imbalance)
# library(smotefamily)


# input data --------------------------------------------------------------

# clinical data
outcome = 'died'

getwd()
setwd("/mnt/std-pool/homedirs/lzhang/coviRNA/")
clinic <- read.xlsx('coviRNA/data/clinic_outcomes_imp_splitDied.xlsx')
dim(clinic)
summary(as.factor(clinic$cohort))
summary(as.factor(clinic[,outcome]))

var_clinic <- c('age','sex_male', 'smoker_ex', 'smoker_curr')
summary(clinic[,var_clinic])


# clinical subset (cohorts) and with outcome
cohorts <- c('NAPKON','LIH','UEDIN')

clinic_sub <- clinic[clinic$cohort %in% cohorts,]
dim(clinic_sub)
summary(as.factor(clinic_sub[,outcome]))
rownames(clinic_sub) <- clinic_sub$firalisID
clinic_sub$group <- as.factor(clinic_sub[,outcome])
summary(clinic_sub$group)

# row counts
load('coviRNA/counts_paxgene.RData')
ls()

dim(counts_all)
counts_all[1:5, 1:5]

counts_sub <- counts_all[fimics_new, clinic_sub$firalisID]
dim(counts_sub)

counts_sub <- counts_sub[rowSums(counts_sub >= 5) >= 3,]
dim(counts_sub)

counts_sub$seqID <- rownames(counts_sub)

# voom transformation
library(limma)

counts_3cohorts <- counts_all[,colnames(counts_all) %in% napkon_uedin_lih$firalisID]
dim(counts_3cohorts)
counts_3cohorts <- counts_3cohorts[fimics_new,]
dim(counts_3cohorts)
counts_3cohorts <- counts_3cohorts[rowSums(counts_3cohorts >= 5) >= 3,]
dim(counts_3cohorts)

v <- voom(counts_sub)
df_v <- as.data.frame(v$E)
dim(df_v)
df_v$seqID <- rownames(df_v)
getwd()
write.csv(df_v, file = 'voom_napkon_556_reads5samples3.csv', row.names = F, quote = F)

# voom transformed data
voom <- read.csv('coviRNA/data/voom_napkon_uedin_lih_reads5samples3.csv')

dim(voom)
voom[1:5,1:5]
rownames(voom) <- voom$seqID
lnc_all <- voom$seqID
voom$seqID <- NULL

voom_t <- as.data.frame(t(voom))
voom_t[1:5,1:5]
voom_t$firalisID <- rownames(voom_t)

merged <- merge(clinic_sub, voom_t, by = 'firalisID')
dim(merged)
rownames(merged) <- merged$firalisID
merged$group <- as.factor(merged[,outcome])
summary(merged$group)
summary(as.numeric(merged$smoker_curr))

merged$sex_male <- as.factor(merged$sex_male)
merged$smoker_ex <- as.factor(merged$smoker_ex)
merged$smoker_curr <- as.factor(merged$smoker_curr)

dim(merged)

# sampling for the same number of alive patients
r <- 5 ##### threshold of reads

samples_1_ <- clinic_sub[clinic_sub[,outcome] == 1, "firalisID"]
samples_0_ <- clinic_sub[clinic_sub[,outcome] == 0,"firalisID"]
length(samples_1_)
length(samples_0_)

test_1 <- round(62*0.2, digits = 0)
test_1
train_1 <- 62-test_1
train_1

lst_lnc <- list()

# when DESeq2 has been down
lst_lnc <- readRDS('deseq2_napkon_uedin_lih_died_08_100.rds')
length(lst_lnc)
class(lst_lnc[[1]])
head(lst_lnc[[1]])

# output lists
lst_boruta_confirmed <- list()
lst_boruta_conf_tent <- list()

lst_vif_confirmed <- list()
lst_vif_conf_tent <- list()

# lst_wilcox <- list()

df_resample <- as.data.frame(matrix(nrow = nrow(clinic_sub), ncol = 101,data = 2))
colnames(df_resample) <- c('firalisID', paste0('train_', seq(1,100,1)))
colnames(df_resample)
df_resample$firalisID <- clinic_sub$firalisID
df_resample[1:5,1:5]
dim(df_resample)

# resampled index
df_resample <- read.csv('coviRNA/resample_balanced_napkon_uedin_lih_08_100.csv')
dim(df_resample)
summary(as.factor(df_resample$train_1))
summary(as.factor(clinic_sub[df_resample[df_resample$train_1 == 1,"firalisID"],"died"]))


# v <- 'p'
padj <- 0.00001


for (i in c(1:100)) {
  time_start <- Sys.time()
  print(paste0('#####################',i,'###################'))
  
  # ##### data split #####
  # # set seeds for reproducibility
  # set.seed(i)
  # # samples_1_ <- dat_use[dat_use[,outcome] == 1,"firalisID"]
  # samples_1 <- sample(samples_1_, size = train_1, replace = F) # training group 1
  # # samples_1_test <- samples_1_[!(samples_1_ %in% samples_1)]
  # # samples_1 <- samples_1_
  # print(length(samples_1))
  # 
  # # samples_0_ <- dat_use[dat_use[,outcome] == 0,"firalisID"]
  # samples_0 <- sample(samples_0_, size = train_1, replace = F) # training group 0
  # # samples_0_test <- sample(samples_0_[!(samples_0_ %in% samples_0)], size = length(samples_1_test))
  # print(length(samples_0))
  
  # from outside
  samples_train <- df_resample[df_resample[,paste0('train_',i)] == 1, "firalisID"]
  # samples_0 <- df_resample[df_resample[,paste0('train_',i)] == 0, "firalisID"]

  # length(samples_1)
  # length(samples_0)
  # length(samples_1_test)
  # length(samples_0_test)
  
  # df_resample[df_resample$firalisID %in% c(samples_0, samples_1), paste0('train_',i)] <- 1
  # df_resample[df_resample$firalisID %in% c(samples_0_test, samples_1_test), paste0('train_',i)] <- 0
  #####
  
  #### DE lncRNAs #####
  # # DESeq2
  dat_counts <- counts_sub[, samples_train]
  dat_cond <- clinic_sub[samples_train,]
  
  samples_0 <- dat_cond[dat_cond$group == 0,"firalisID"]
  samples_1 <- dat_cond[dat_cond$group == 1,"firalisID"]
  #
  # # filter low expressed lncRNAs
  is.1 <- rowSums(dat_counts[,samples_0] >= r) >= length(samples_0)/2
  is.2 <- rowSums(dat_counts[,samples_1] >= r) >= length(samples_1)/2
  # lnc_all <- rownames(dat_counts[is.1 | is.2,])
  # length(lnc_all)
  dat_counts_ <- dat_counts[is.1 | is.2,]
  print(dim(dat_counts_))

  dds_filter <- DESeqDataSetFromMatrix(countData = dat_counts_,
                                       colData = dat_cond,
                                       design = ~group)
  dds_filter <- DESeq(dds_filter, test = 'Wald')
  res_filter <- results(dds_filter,
                        contrast = c('group',"1","0"),
                        alpha = 0.05)
  lst_lnc[[i]] <- res_filter
  lnc_de <- rownames(res_filter[which(res_filter$padj < padj),]) #####
  # 
  # # when DESeq2 has been done
  # df_lnc <- lst_lnc[[i]]
  # lnc_de <- rownames(df_lnc[which(df_lnc$padj < padj),])
  print(length(lnc_de))
  # 
  #####
  
  # # merged_$firalisID <- rownames(merged_)
  dat_merged <- merged[samples_train, c(var_clinic,lnc_de, 'group')]
  print(dim(dat_merged))
  dat_merged[,c('age',lnc_de)] <- scale(dat_merged[,c('age',lnc_de)])
  
  ##### boruta feature selection #####
  
  set.seed(i)
  res_boruta <- Boruta(x = dat_merged[,c(var_clinic,lnc_de)], y = dat_merged$group, doTrace = 0, maxRuns = 500)
  var_boruta_2 <- names(res_boruta$finalDecision[res_boruta$finalDecision %in%  c('Confirmed','Tentative')]) #####
  var_boruta <- names(res_boruta$finalDecision[res_boruta$finalDecision == 'Confirmed']) #####
  print(length(var_boruta))
  lst_boruta_conf_tent[[i]] <- var_boruta_2
  lst_boruta_confirmed[[i]] <- var_boruta
  #####
   
  ##### multicollinearity #####
  form <- as.formula(paste('group', paste(var_boruta, collapse = ' + '), sep = ' ~ '))
  model <- glm(form, dat_merged, family = 'binomial')
  vif_ <- vif(model)
  var_vif <- names(vif_[vif_ <= 5])
  print(length(var_vif))
  lst_vif_confirmed[[i]] <- var_vif
  
  form <- as.formula(paste('group', paste(var_boruta_2, collapse = ' + '), sep = ' ~ '))
  model <- glm(form, dat_merged, family = 'binomial')
  vif_ <- vif(model)
  var_vif <- names(vif_[vif_ <= 5])
  lst_vif_conf_tent[[i]] <- var_vif
  
  time_end <- Sys.time()
  print(time_end - time_start)
}

length(lst_vif_confirmed)
# lst_vif

# selected times
lst_vif <- lst_vif_confirmed

df_sel_n <- as.data.frame(lst_vif[[1]])
colnames(df_sel_n) <- 'feature'
df_sel_n$n <- rep(1,nrow(df_sel_n))
df_sel_n

for (i in c(2:100)) {
  for (j in lst_vif[[i]]) {
    if (j %in% df_sel_n$feature) {
      df_sel_n[df_sel_n$feature == j,"n"] <- df_sel_n[df_sel_n$feature == j,"n"] + 1
    } else {
      df_sel_n <- rbind(df_sel_n, data.frame(feature=j, n=1))
    }
  }
}

dim(df_sel_n)
View(df_sel_n)

getwd()
saveRDS(lst_vif_confirmed, 'vif10Conf_df000001_napkon_uedin_lih_died_08_100_.rds')
write.xlsx(df_sel_n, 'features_sel_numbers_df000001_napkon_uedin_lih_died_08_100_.xlsx')
write.csv(df_resample, 'resample_balanced_napkon_uedin_lih_batchCorr_died_08_100_.csv', row.names = F, quote = F)

##### data frame selected features at each iteration #####
getwd()
lst_vif <- readRDS('coviRNA/vif5Conf_df000001_napkon_uedin_lih_died_08_100_.rds')

features_1 <- NULL
for (i in c(1:100)) {
  features_ <- lst_vif[[i]]
  features_1 <- unique(c(features_1, features_))
}
length(features_1)

df_features <- as.data.frame(matrix(nrow = 228, ncol = 101))
colnames(df_features) <- c('feature',paste0('select_',seq(1,100,1)))
colnames(df_features)
df_features$feature <- features_1
for (i in c(1:100)) {
  features_ <- lst_vif[[i]]
  df_features[,paste0('select_',i)] <- 0
  df_features[df_features$feature %in% features_, paste0('select_',i)] <- 1
}
dim(df_features)

write.csv(df_features, 'features_100iter_lncDE000001.csv', row.names = F, quote = F)

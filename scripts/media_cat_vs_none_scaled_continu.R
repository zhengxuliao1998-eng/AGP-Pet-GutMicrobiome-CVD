# devtools::install_github("mbsohn/cmmb")
library(cmmb)
library(phyloseq)
library(tidyverse)


meta <- read.table("/N/project/Watson_AGP/Zhengxu/DataSet/cleaned/meta/meta_07202026.txt")
load("/N/u/liaozh/Quartz/Desktop/Link to Zhengxu/DataSet/cleaned/biom/biom_07202026.Rdata")
biom@otu_table <- otu_table(subset(biom@otu_table,select = meta$ID))

# select 0.1 taxon
ind <- which(apply(biom@otu_table,1,function(x){
  mean(x>0)
})>=0.1)
biom@otu_table <- biom@otu_table[ind,]

taxa_id_set <- rownames(as.data.frame(biom@tax_table[which(biom@tax_table[,3]%in%c("c__Mollicutes", "c__Actinobacteria")),]))
biom@otu_table <- biom@otu_table[rownames(biom@otu_table)%in%taxa_id_set,]

# select meta info
meta_design <- meta%>%
  select(cat_only, dog_only, cat_and_dog,
         age, 
         age_group55, 
         sleep_duration_regroup,
         sex ,
         race_recode ,
         smoker , 
         exercise_recode ,
         alcohol_recode ,
         education_recode ,
         bmi , 
         bmi_cat,
         protein_intake ,
         antibiotic_recent_take
  )%>%
  mutate(age_group55 = as.factor(age_group55),
         age_group55 = relevel(age_group55, ref = "under55"),
         sleep_duration_regroup = as.factor(sleep_duration_regroup),
         sex = as.factor(sex),
         race_recode = as.factor(race_recode),
         education_recode = as.factor(education_recode),
         protein_intake = as.factor(protein_intake),
         antibiotic_recent_take= as.factor(antibiotic_recent_take)
  )

# select meta info
meta_design <- data.frame(
  model.matrix(
    ~cat_only+dog_only+cat_and_dog+
      age + 
      sleep_duration_regroup +
      sex +
      race_recode +
      smoker + 
      exercise_recode +
      alcohol_recode +
      education_recode +
      protein_intake +
      antibiotic_recent_take
    ,
    meta_design))[,-1]

meta_design <- meta_design%>%mutate(across(-c(cat_only,dog_only,cat_and_dog), scale))
meta_design <- as.data.frame(meta_design)
otu <- biom@otu_table
otu[which(otu==0)] = 0.5
remove(biom)

# separate set and save to different data.frame
### Part 1: separate otu tables
cat_only_otu = t(subset(otu,select = meta$ID[meta$cat_only==1]))
# dog_only_otu = t(subset(otu,select = meta$ID[meta$dog_only==1]))
# cat_dog_otu = t(subset(otu,select = meta$ID[meta$cat_and_dog==1]))
no_pet_otu = t(subset(otu,select = meta$ID[meta$cat_only==0&meta$dog_only==0&meta$cat_and_dog==0]))

### Part 2: separate meta dtables
cat_only_meta = meta_design%>%
  filter(cat_only==1)%>%select(-cat_only,-dog_only,-cat_and_dog)
# dog_only_meta = meta_design%>%
#   filter(dog_only==1)%>%select(-cat_only,-dog_only,-cat_and_dog)
# cat_dog_meta = meta_design%>%
#   filter(cat_and_dog==1)%>%select(-cat_only,-dog_only,-cat_and_dog)
no_pet_meta = meta_design%>%
  filter(cat_only==0&dog_only==0&cat_and_dog==0)%>%select(-cat_only,-dog_only,-cat_and_dog)

### Part 3: Separate outcome table
cat_only_Y = meta%>%filter(cat_only==1)%>%select(cvd)
# dog_only_Y = meta%>%filter(dog_only==1)%>%select(cvd)
# cat_dog_Y = meta%>%filter(cat_and_dog==1)%>%select(cvd)
no_pet_Y = meta%>%filter(cat_only==0&dog_only==0&cat_and_dog==0)%>%select(cvd)

### Part 4: Separate treatment vector
cat_only_tr = meta%>%filter(cat_only==1)%>%select(cat_only,dog_only,cat_and_dog,no_pet)
# dog_only_tr = meta%>%filter(dog_only==1)%>%select(cat_only,dog_only,cat_and_dog,no_pet)
# cat_dog_tr = meta%>%filter(cat_and_dog==1)%>%select(cat_only,dog_only,cat_and_dog,no_pet)
no_pet_tr = meta%>%filter(cat_only==0&dog_only==0&cat_and_dog==0)%>%select(cat_only,dog_only,cat_and_dog,no_pet)


cat_only_otu <- cat_only_otu/apply(cat_only_otu,1,sum)
# dog_only_otu <- dog_only_otu/apply(dog_only_otu,1,sum)
# cat_dog_otu <- cat_dog_otu/apply(cat_dog_otu,1,sum)
no_pet_otu <- no_pet_otu/apply(no_pet_otu,1,sum)

rownames(cat_only_otu) <- NULL
# rownames(dog_only_otu) <- NULL
# rownames(cat_dog_otu) <- NULL
rownames(no_pet_otu) <- NULL

# characters are subject to be factors
# note: no 0 in otu compositional data.
Y=c(cat_only_Y$cvd,no_pet_Y$cvd)
M=rbind(as.matrix(cat_only_otu),as.matrix(no_pet_otu))
tr=c(cat_only_tr$cat_only,no_pet_tr$cat_only)
X=as.matrix(sapply(rbind(cat_only_meta,no_pet_meta), as.numeric))

write.csv(M,"/N/u/liaozh/Quartz/Desktop/Link to Zhengxu/Analysis/Rprogram/job_scripts/temp_files/cat_vs_none/M.txt")
write.csv(X,"/N/u/liaozh/Quartz/Desktop/Link to Zhengxu/Analysis/Rprogram/job_scripts/temp_files/cat_vs_none/X.txt")
write.csv(Y,"/N/u/liaozh/Quartz/Desktop/Link to Zhengxu/Analysis/Rprogram/job_scripts/temp_files/cat_vs_none/Y.txt")
write.csv(tr,"/N/u/liaozh/Quartz/Desktop/Link to Zhengxu/Analysis/Rprogram/job_scripts/temp_files/cat_vs_none/tr.txt")

M <- as.matrix(read.csv("M.txt")[,-1])
X <- as.matrix(read.csv("X.txt")[,-1])
Y <- read.csv("Y.txt")[,-1]
tr <- read.csv("tr.txt")[,-1]

st <- Sys.time()
cat_vs_none_media_result <- cmmb(Y=Y,
                                 M=M,
                                 tr=tr,
                                 X=X,
                                 p.value = T, n.cores = 16, n.boot = 1000) 

save(cat_vs_none_media_result, file = "cat_vs_none_bo.Rdata")
tm.cost <- Sys.time()-st
tm.cost

# Pre-requisite
#install.packages("magrittr")
#install.packages("ggolot2")
#nstall.packages("ggridges")
#install.packages("readxl")

# Author: Peter Lin
# Last update: 2024.11.19

# 0. Environment------
library(magrittr)
library(ggplot2)
library(ggridges)
setwd("C:/Users/User/Desktop/R1152/Tomato RADseq x ML/accession mapping/full_datasets")

# 1. Load data--------
Sys.setlocale("LC_ALL", "Chinese (Traditional)_Taiwan.950")
#tb_2018_cool <- read.csv("107種源番茄 (19號).csv", sep = ",",
#                           fileEncoding = "UTF16")
#read.csv() resulted partial trucated input for unknown reason
tb <- readxl::read_xlsx("2018_cool.xlsx",
                        col_types = c("date", rep("text",3), "numeric", "skip",
                                      "numeric", "text", rep("numeric",3), "text", rep("numeric",2), "text",
                                      rep("numeric",15), "skip", "skip"))
# Skipped columns are "T4.處理", "X1.備註(1)", "X2.備註(2)"
# These columns are empty

# Correct mis-typed decimal points
extraDeci <- grep("6..35", tb$`P5.糖度(brix)`)
tb$`P5.糖度(brix)`[extraDeci] <- 6.35
extraDeci <- grep("8..22", tb$`P5.糖度(brix)`)
tb$`P5.糖度(brix)`[extraDeci] <- 8.22

# Delete ambiguous records
ambi <- grep("[(]3.8[)]", tb$`P5.糖度(brix)`)
tb$`P5.糖度(brix)`[ambi] <- 0
ambi <- grep("[(]4[)]", tb$`P5.糖度(brix)`)
tb$`P5.糖度(brix)`[ambi] <- 0
ambi <- grep("8[?]", tb$`P5.糖度(brix)`)
tb$`P5.糖度(brix)`[ambi] <- 0

# Convert brix record into numeric
tb$`P5.糖度(brix)` <- as.numeric(tb$`P5.糖度(brix)`)

map <- readxl::read_xlsx("accession mapping.xlsx")
comp_phen <- map$核心種原編號[map$COMP_PHEN == 1] #有完整平均值資料的品系
big_F <- comp_phen %in% map$核心種原編號[map$BIG_F == 1]
small_F <- comp_phen %in% map$核心種原編號[map$BIG_F == 0]

# 2. filtering----------
#Filter out unwanted accessions
acc_num_len <- nchar(tb$T2.品系代碼)
table(acc_num_len)

keep <- acc_num_len <= 4
#康樂信件敘述:
#核心種原的品種代碼都是三位數以內的
#如果後面有A或B表示不同重複
tb_filtered <- tb[keep,]
tb_filtered$T2.品系代碼 <- as.numeric(tb_filtered$T2.品系代碼)

unique(tb_filtered$T2.品系代碼) %>% sort() #289 accessions

tb_filtered <- tb_filtered[tb_filtered$T2.品系代碼 %in% comp_phen,]
length(unique(tb_filtered$T2.品系代碼)) # 260 accessions


# screen for suspicious missing data or zeros
any.zero <- function(x){
  sum(x == 0, na.rm = F)
}
#out <- apply(tb_filtered[,c(8:13,15:17)], MARGIN = 1, FUN = any.zero)
out <- apply(tb_filtered[,c(8:9)], MARGIN = 1, FUN = any.zero)
out <- out > 0 | is.na(out)
tb_out <- tb_filtered[out,] #fruits with one or more missing values
tb_filtered <- tb_filtered[out == F,c(1:11,14:17,13)] #只留下 長度、寬度、果重、糖度、果色LAB、體積
#writexl::write_xlsx(tb_out, path = "有疑問的資料_107.xlsx")
#tb_filtered <- tb_filtered[out != T,]

#drop records with fruit length/width < 0.2
drop <- tb_filtered$`P2.長度(cm)` < 0.2 | tb_filtered$`P3.寬度(cm)` < 0.2
sum(drop, na.rm = T) # 21
tb_filtered <- tb_filtered[drop == F,]

#assign NA for every missing measurement
for (i in c(8:11,13:15)){
  NAs <- c()
  NAs <- tb_filtered[,i] == 0 | is.na(tb_filtered[,i])
  NAs <- NAs[,1]
  tb_filtered[NAs,i] <- NA
}

table(tb_filtered$`C1.整體顏色(原果色)`)
#咖啡       紅     粉紅     淡黃       紫       黃     黃紅       綠   綠帶紫 綠帶紫黑 
#48    40046     3041      774      557     2219       75     2760     2286       99 
#綠紫       橘     橘紅 
#583      834        8 


#calculate the lycopene content with a/b
tb_filtered$lyco_1 <- 11.848*(tb_filtered$`C3.顏色-A`/tb_filtered$`C4.顏色-B`) + 1.5471
tb_filtered$lyco_1[tb_filtered$lyco_1 > 30 | tb_filtered$lyco_1 < 0] <- NA

#calculate the lycopene content with (a/b)^2 (not used)
#tb_filtered$lyco_2 <- 8.7073*c((tb_filtered$`C3.顏色-A`/tb_filtered$`C4.顏色-B`)^2) + 1.5212
#tb_filtered$lyco_2[tb_filtered$lyco_2 > 50 ] <- NA


# Fruits WITH effective lycopene content estimation (count by fruit color)
tb_filtered$`C1.整體顏色(原果色)`[is.na(tb_filtered$lyco_1)==F] %>% table()
#咖啡       紅     粉紅     淡黃       紫       黃     黃紅      綠   綠帶紫 綠帶紫黑 
#45    35915     2747      754       98     2213       75     1870     1708       96 
#綠紫       橘     橘紅 
#348      834        8

# Fruits WITHOUT effective lycopene content estimation (count by fruit color)
tb_filtered$`C1.整體顏色(原果色)`[is.na(tb_filtered$lyco_1)] %>% table()
#咖啡       紅     粉紅     淡黃       紫       黃       綠   綠帶紫 綠帶紫黑     綠紫 
#3     4131      294       20      459        6      890      578        3      235 

nolyco <- is.na(tb_filtered$lyco_1)
red <- tb_filtered$`C1.整體顏色(原果色)`=="紅"
pink <- tb_filtered$`C1.整體顏色(原果色)`=="粉紅"
yellow <- tb_filtered$`C1.整體顏色(原果色)`=="黃"
lyellow <- tb_filtered$`C1.整體顏色(原果色)`=="淡黃"
green <- c(1:nrow(tb_filtered)) %in% grep("綠", tb_filtered$`C1.整體顏色(原果色)`)
purple <- tb_filtered$`C1.整體顏色(原果色)`=="紫"
brown <- tb_filtered$`C1.整體顏色(原果色)`=="咖啡"
is.na(tb_filtered$`C1.整體顏色(原果色)`) %>% sum() #47
withcol <- is.na(tb_filtered$`C1.整體顏色(原果色)`) == F

#Fix ineffective lycopene values for "紅" fruits with median
summary(tb_filtered$lyco_1[red & withcol])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.129  18.198  21.062  20.887  23.837  29.998    4131 
tb_filtered$lyco_1[red & withcol & nolyco] <- mean(tb_filtered$lyco_1[red & withcol], na.rm = T)

#Fix ineffective lycopene values for "粉紅" fruits with median
summary(tb_filtered$lyco_1[pink & withcol])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#3.865  20.502  23.647  23.020  26.211  29.977     294 
tb_filtered$lyco_1[pink & withcol & nolyco] <- mean(tb_filtered$lyco_1[pink & withcol], na.rm = T)

#Fix ineffective lycopene values for "黃" fruits with median
summary(tb_filtered$lyco_1[yellow & withcol])
#Min.   1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.0288  4.6656  7.0860  6.6757  8.4614 13.0982       6
tb_filtered$lyco_1[yellow & withcol & nolyco] <- mean(tb_filtered$lyco_1[yellow & withcol], na.rm = T)

#Fix ineffective lycopene values for "淡黃" fruits with median
summary(tb_filtered$lyco_1[lyellow & withcol])
#Min.     1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#0.001862 2.154698 3.111960 3.091842 3.885739 7.932522       20 
tb_filtered$lyco_1[lyellow & withcol & nolyco] <- mean(tb_filtered$lyco_1[lyellow & withcol], na.rm = T)

#Fix ineffective lycopene values for "綠" fruits with median
summary(tb_filtered$lyco_1[green & withcol])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0017  1.3377  3.4370  5.5726  7.6215 29.9093    1706
tb_filtered$lyco_1[green & withcol & nolyco] <- mean(tb_filtered$lyco_1[green & withcol], na.rm = T)

#tb_purple <- tb_filtered[brown,]

#分離果色含"綠"的紀錄
#GandP_tb <- tb_filtered[grep("綠",tb_filtered$`C1.整體顏色(原果色)`),]
#GandP_tb <- GandP_tb[c(is.na(GandP_tb$lyco_1) | is.na(GandP_tb$lyco_2)) == F,]
#writexl::write_xlsx(GandP_tb, path = "green & purple records.xlsx")

# plots of estimated lycopene------------
lyco_tb <- tb_filtered[,12:17]
lyco_tb <- lyco_tb[is.na(lyco_tb$lyco_1) == F,]
unique(lyco_tb$`C1.整體顏色(原果色)`) %>% sort()

# Visualize the distributions of lycopene content of each fruit color
ggplot(data = lyco_tb)+
  geom_density_ridges(mapping = aes(x = lyco_1, y = `C1.整體顏色(原果色)`,
                                    fill = `C1.整體顏色(原果色)`), alpha = 0.5)+
  scale_fill_manual(breaks = c("綠",      
                               "綠帶紫","綠帶紫黑","綠紫","紫",
                               "咖啡", "紅","粉紅","橘紅","黃紅",
                               "橘","黃","淡黃"),
                    values = c("#6fb510",
                               "#087d46","#3b224f","#087d46","#a763dc",
                               "#954e07","#f1280c","#f2a6a6","#df520b","#f1ad71",
                               "#f2b52a","#fbd82a","#fdf3c1"))+
  xlab(label = "Lycopene content by a/b (mg/100g)")+
  xlim(c(0,50))
  
#ggplot(data = lyco_tb)+
#  geom_point(mapping = aes(x = lyco_1, y = lyco_2, color = `C1.整體顏色(原果色)`))+
#  scale_color_manual(breaks = c("綠",      
#                                "綠帶紫","綠帶紫黑","綠紫","紫",
#                                "咖啡", "紅","粉紅","橘紅","黃紅",
#                                "橘","黃","淡黃"),
#                     values = c("#6fb510",
#                                "#087d46","#3b224f","#087d46","#a763dc",
#                                 "#954e07","#f1280c","#f2a6a6","#df520b","#f1ad71",
#                                 "#f2b52a","#fbd82a","#fdf3c1"))+
#   xlab(label = "Lycopene content by a/b (mg/100g)")+
#   ylab(label = "Lycopene content by (a/b)^2 (mg/100g)")
# 
# ggplot(data = lyco_tb)+
#   geom_density_ridges(mapping = aes(x = lyco_2, y = `C1.整體顏色(原果色)`,
#                                     fill = `C1.整體顏色(原果色)`), alpha = 0.5)+
#   scale_fill_manual(breaks = c("綠",      
#                                "綠帶紫","綠帶紫黑","綠紫","紫",
#                                "咖啡", "紅","粉紅","橘紅","黃紅",
#                                "橘","黃","淡黃"),
#                     values = c("#6fb510",
#                                "#087d46","#3b224f","#087d46","#a763dc",
#                                "#954e07","#f1280c","#f2a6a6","#df520b","#f1ad71",
#                                "#f2b52a","#fbd82a","#fdf3c1"))+
#   xlab(label = "Lycopene content by (a/b)^2 (mg/100g)")+
#   xlim(c(0,50))
#   
# ggplot(data = lyco_tb)+
#   geom_point(mapping = aes(x = lyco_1, y = lyco_2,
#                            color = `C1.整體顏色(原果色)`, alpha = `C1.整體顏色(原果色)`))+
#   scale_color_manual(breaks = c("綠",      
#                                 "綠帶紫","綠帶紫黑","綠紫","紫",
#                                 "咖啡", "紅","粉紅","橘紅","黃紅",
#                                 "橘","黃","淡黃"),
#                      values = c("#6fb510",
#                                 "#087d46","#3b224f","#087d46","#a763dc",
#                                 "#954e07","#f1280c","#f2a6a6","#df520b","#f1ad71",
#                                 "#f2b52a","#fbd82a","#fdf3c1"))+
#   scale_alpha_manual(breaks = c("綠",      
#                                 "綠帶紫","綠帶紫黑","綠紫","紫",
#                                 "咖啡", "紅","粉紅","橘紅","黃紅",
#                                 "橘","黃","淡黃"),
#                      values = c(1,
#                                 1,1,1,1,
#                                 0,0,0,0,0,
#                                 0,0,0))+
#   xlab(label = "Lycopene content by a/b (mg/100g)")+
#   ylab(label = "Lycopene content by (a/b)^2 (mg/100g)")
# 
# ggplot(data = lyco_tb)+
#   geom_point(mapping = aes(x = lyco_1, y = lyco_2,
#                            color = `C1.整體顏色(原果色)`, alpha = `C1.整體顏色(原果色)`))+
#   scale_color_manual(breaks = c("綠",      
#                                 "綠帶紫","綠帶紫黑","綠紫","紫",
#                                 "咖啡", "紅","粉紅","橘紅","黃紅",
#                                 "橘","黃","淡黃"),
#                      values = c("#6fb510",
#                                 "#087d46","#3b224f","#087d46","#a763dc",
#                                 "#954e07","#f1280c","#f2a6a6","#df520b","#f1ad71",
#                                 "#f2b52a","#fbd82a","#fdf3c1"))+
#   scale_alpha_manual(breaks = c("綠",      
#                                 "綠帶紫","綠帶紫黑","綠紫","紫",
#                                 "咖啡", "紅","粉紅","橘紅","黃紅",
#                                 "橘","黃","淡黃"),
#                      values = c(0,
#                                 0,0,0,0,
#                                 1,1,1,1,1,
#                                 1,1,1))+
#   xlab(label = "Lycopene content by a/b (mg/100g)")+
#   ylab(label = "Lycopene content by (a/b)^2 (mg/100g)")

# Screen for unusual fruit weight---------
WVratio <- tb_filtered$`P4.果重(g)`/tb_filtered$`P7.體積(cm3)` #Weight-Volume ratio
sum(WVratio >= 1.2, na.rm = T) #9265
sum(c(WVratio < 1.2 & WVratio > 1.1), na.rm = T) #3182
sum(c(WVratio >= 0.9 & WVratio <= 1.1), na.rm = T) #28223
sum(c(WVratio < 0.9 & WVratio > 0.8), na.rm = T) #9244
sum(WVratio <= 0.8, na.rm = T) #2735
summary(WVratio) #Median = 0.9852
med <- 0.9852

sum(WVratio >= med+0.2, na.rm = T) #9568
sum(c(WVratio < med+0.2 & WVratio > med+0.1), na.rm = T) #3811
sum(c(WVratio >= med-0.1 & WVratio <= med+0.1), na.rm = T) #29151
sum(c(WVratio < med-0.1 & WVratio > med-0.2), na.rm = T) #7933
sum(WVratio <= med-0.2, na.rm = T) #2186

ggplot(data=data.frame(WVratio = WVratio))+
  geom_density(mapping = aes(x=WVratio))+
  xlim(c(0,2))+
  xlab("Weight-Volume ratio")

correctW <- tb_filtered$`P4.果重(g)`
toCorrect <- WVratio > med+0.2
toCorrect[is.na(toCorrect)] <- FALSE
correctW[toCorrect] <- tb_filtered$`P7.體積(cm3)`[toCorrect]*(med+0.2)
toCorrect <- WVratio < med-0.2
toCorrect[is.na(toCorrect)] <- FALSE
correctW[toCorrect] <- tb_filtered$`P7.體積(cm3)`[toCorrect]*(med-0.2)
correctW[is.na(WVratio)] <- tb_filtered$`P7.體積(cm3)`[is.na(WVratio)]*med
summary(correctW)

ggplot(data=data.frame(Weight=c(tb_filtered$`P4.果重(g)`,correctW),
                       Group=rep(c("Raw","Corrected"),each=length(correctW))))+
  geom_density(mapping = aes(x=Weight, fill=Group), alpha=0.4)+
  xlim(c(0,50))+
  xlab("Fruit weight(g)")

# Calculate fruit weight per plant (fwpp)--------
a <- c() # a for accession
p <- c() # p for plant count
fwpp <- c() # fwpp for fruit weight per plant
pass <- c() #check if every plant has 3 or more fruits
for (i in comp_phen){
  row <- tb_filtered$T2.品系代碼 == i #obtain fruits belong to accession i
  dt <- tb_filtered[row,] #all fruits belong to accession i
  n <- unique(dt$T3.株號) %>% length() #number of plants of accession i
  a <- c(a,i)
  p <- c(p,n)
  av <- table(dt$T3.株號) >= 3 #check if every plant has 3 or more fruits
  av <- sum(av) == n #check if every plant has 3 or more fruits
  pass <- c(pass, av)
  fw <- sum(dt$`P4.果重(g)`, na.rm = T) #weight of all fruits of accession i
  fw <- fw / length(unique(dt$T3.株號[is.na(dt$`P4.果重(g)`)==F]))
        #weight of all fruits of accession i / #number of plants with at least one "weighted" fruit
  if (length(unique(dt$T3.株號[is.na(dt$`P4.果重(g)`)==F])) > 6){
    fw <- fw*length(unique(dt$T3.株號[is.na(dt$`P4.果重(g)`)==F]))/6
  }
  fwpp <- c(fwpp, fw)
}
pc_tb <- data.frame(accession = a, plant_count = p, FWPP = fwpp)

#calculate means, sds, cvs for traits within each population
cv <- c()
avg <- c()
sl <- c()
for (i in comp_phen){
  row <- tb_filtered$T2.品系代碼 == i
  dt <- tb_filtered[row,]
  sd <- apply(dt[,c(8:11,16)], MARGIN = 2, FUN = sd, na.rm=T)
  mean <- apply(dt[,c(8:11,16)], MARGIN = 2, FUN = mean, na.rm=T)
  sl <- c(sl, mean(unique(dt$`P5.糖度(brix)`), na.rm=T))
  cv <- c(cv, unname(round(sd/mean*100, digits = 2)))
  avg <- c(avg, round(mean,digits = 2))
}

#make a table of coefficient of variations
dim(cv) <- c(5,260)
cv <- data.frame(cv)
colnames(cv) <- comp_phen
rownames(cv) <- c("fruit length", "fruit width", "fruit weight",
                  "sugar level", "estimated lycopene")

#view the distribution of coefficient of variations across accessions
plot(1:ncol(cv),sort(cv[3,]), ylim = c(0,100)) # Fruit weight
plot(1:ncol(cv),sort(cv[5,]), ylim = c(0,100)) # Lycopene content
hist(as.numeric(cv[3,]), main = "Coefficient of Variation of Fruit Weight \nwithin Each Accession (2018 cool)",
     ylab = "Number of Accessions", xlab = "Coefficient of Variation (%)",
     xlim = c(0,max(cv[3,])), breaks = 20)

hist(as.numeric(cv[4,]), main = "Coefficient of Variation of Sugar Level \nwithin Each Accession (2018 cool)",
     ylab = "Number of Accessions", xlab = "Coefficient of Variation (%)",
     xlim = c(0,max(cv[4,])), breaks = 20)

#make a table of trait averages
dim(avg) <- c(5, 260)
avg <- data.frame(avg)
colnames(avg) <- comp_phen
rownames(avg) <- c("fruit length", "fruit width", "fruit weight",
                  "sugar level","estimated lycopene")
avg <- data.frame(t(avg))
avg$serial_num <- as.numeric(rownames(avg))
avg <- avg[,c(6,1:5)]

avg$sugar.level <- sl

#compare the mean fruit weight provided by TARI and the mean value calculated from their raw data
avg_2018 <- readxl::read_xlsx(path = "../phenotype data.xlsx",
                              sheet = "2018-2019")
avg_2018 <- avg_2018[,c(1,17,18,19,21,20)]
avg_2018 <- avg_2018[avg_2018$serial_num %in% comp_phen,]
#mean difference between fruit weight calculated from raw data and the provided fruit weight
c(avg$fruit.weight - avg_2018$`fruit_weight(g)_2018cool`) %>% mean() #0.6053
#numbers of accessions that have fruit weight calculated from raw data identical with the provided fruit weight
c(avg$fruit.weight == avg_2018$`fruit_weight(g)_2018cool`) %>% sum() #237
unique(tb_out$T2.品系代碼) %>% length() #4 accessions 

#writexl::write_xlsx(avg_2018, path = "avg_2018_fromTARI.xlsx")
#writexl::write_xlsx(avg, path = "avg_2018_fromRAW.xlsx")

#calculate yields and compare them with the yields provided by TARI
yield <- c(pc_tb$FWPP)*33000/1000 #fruit weight per plant*33000/1000
plot(yield, avg_2018$`yield(kg/ha)_2018cool`)
cor(yield, avg_2018$`yield(kg/ha)_2018cool`) #0.9117624
sum(yield == avg_2018$`yield(kg/ha)_2018cool`)
yield_comp <- data.frame(accession = comp_phen,
                         yield_from_RAW = yield,
                         yield_2018cool_TARI = avg_2018$`yield(kg/ha)_2018cool`,
                         diff = yield - avg_2018$`yield(kg/ha)_2018cool`)
yield_comp$diff_percent <- yield_comp$diff/yield_comp$yield_2018cool_TARI*100

#writexl::write_xlsx(yield_comp, path = "產量比較_107 涼季.xlsx")

avg$yield <- yield_comp$yield_from_RAW

avg_fromRAW_2018 <- avg

save(avg_fromRAW_2018, file = "avg_fromRAW_2018.RData")

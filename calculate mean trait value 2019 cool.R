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


# 1. load data--------
Sys.setlocale("LC_ALL", "Chinese (Traditional)_Taiwan.950")
#tb_2018_cool <- read.csv("107種源番茄 (19號).csv", sep = ",",
#                           fileEncoding = "UTF16")
#read.csv() resulted partial trucated input for unknown reason
tb <- readxl::read_xlsx("2019_fall.xlsx",
                        col_types = c(rep("guess",5), "skip",
                                      rep("guess", 5), "text", #Column 12 "P5.糖度(brix)" contains description and numbers
                                      rep("guess", 18), "skip", "skip"))
# Skipped columns are "T4.處理", "X1.備註(1)", "X2.備註(2)"
# These columns are empty

# Delete ambiguous records
unique(tb$`P5.糖度(brix)`) %>% sort(decreasing = T)
# 3 kinds of descriptions were found among other measurements
# "綠果"               "無法測"             "全綠果"
# These descriptions are replaced with NA
NA_brix <- c(grep("綠果", tb$`P5.糖度(brix)`),
             grep("無法測", tb$`P5.糖度(brix)`),
             grep("全綠果", tb$`P5.糖度(brix)`))
tb$`P5.糖度(brix)`[NA_brix] <- NA
tb$`P5.糖度(brix)` <- as.numeric(tb$`P5.糖度(brix)`)

map <- readxl::read_xlsx("accession mapping.xlsx")
comp_phe <- map$核心種原編號[map$COMP_PHEN == 1]
big_F <- comp_phe %in% map$核心種原編號[map$BIG_F == 1]
small_F <- comp_phe %in% map$核心種原編號[map$BIG_F == 0]

# 2. filtering----------
#Filter out unwanted accessions and unify different batches of the same accessions
acc_num_len <- nchar(tb$T2.品系代碼)
table(acc_num_len)
keep <- acc_num_len <= 4
sum(keep) # 104649

tb_filtered <- tb[keep,]
unique(tb_filtered$T2.品系代碼) %>% sort() 
# "150A", "150B" contain English letter to distinguish between replicates
acc150 <- c(grep("150A", tb_filtered$T2.品系代碼), grep("150B", tb_filtered$T2.品系代碼))
tb_filtered$T2.品系代碼[acc150] <- "150"
# "196A", "196B" contain English letter to distinguish between replicates
acc196 <- c(grep("196A", tb_filtered$T2.品系代碼), grep("196B", tb_filtered$T2.品系代碼))
tb_filtered$T2.品系代碼[acc196] <- "196"
# "56A", "56B" contain English letter to distinguish between replicates
acc56 <- c(grep("56A", tb_filtered$T2.品系代碼), grep("56B", tb_filtered$T2.品系代碼))
tb_filtered$T2.品系代碼[acc56] <- "56"
#"81A", "81B" contain English letter to distinguish between replicates
acc81 <- c(grep("81A", tb_filtered$T2.品系代碼), grep("81B", tb_filtered$T2.品系代碼))
tb_filtered$T2.品系代碼[acc81] <- "81"

tb_filtered$T2.品系代碼 <- as.numeric(tb_filtered$T2.品系代碼)

unique(tb_filtered$T2.品系代碼) %>% sort() #288 accessions

tb_filtered <- tb_filtered[tb_filtered$T2.品系代碼 %in% comp_phe,]
length(unique(tb_filtered$T2.品系代碼)) # 260 accessions

# filter out suspicious missing data or zeros
# screen for suspicious missing data or zeros
any.zero <- function(x){
  sum(x == 0, na.rm = F)
}
#out <- apply(tb_filtered[,c(8:13,15:17)], MARGIN = 1, FUN = any.zero)
out <- apply(tb_filtered[,c(8:9)], MARGIN = 1, FUN = any.zero)
out <- out > 0 | is.na(out)
tb_out <- tb_filtered[out,] #fruits with one or more missing values
tb_filtered <- tb_filtered[out==F,c(1:11,14:17,13)] #只留下 長度、寬度、果重、糖度、果色LAB、體積
#writexl::write_xlsx(tb_out, path = "有疑問的資料_108cool.xlsx")
#tb_filtered <- tb_filtered[out != T,]

#drop records with fruit length/width < 0.2
drop <- tb_filtered$`P2.長度(cm)` < 0.2 | tb_filtered$`P3.寬度(cm)` < 0.2
sum(drop, na.rm = T) # 0
tb_filtered <- tb_filtered[drop == F,]

#assign NA for every missing measurement
for (i in c(8:11,13:15)){
  NAs <- c()
  NAs <- tb_filtered[,i] == 0 | is.na(tb_filtered[,i])
  NAs <- NAs[,1]
  tb_filtered[NAs,i] <- NA
}

#calculate the lycopene content with a/b
tb_filtered$lyco_1 <- 11.848*(tb_filtered$`C3.顏色-A`/tb_filtered$`C4.顏色-B`) + 1.5471
tb_filtered$lyco_1[tb_filtered$lyco_1 > 30 | tb_filtered$lyco_1 < 0] <- NA

#calculate the lycopene content with (a/b)^2
#tb_filtered$lyco_2 <- 8.7073*c((tb_filtered$`C3.顏色-A`/tb_filtered$`C4.顏色-B`)^2) + 1.5212
#tb_filtered$lyco_2[tb_filtered$lyco_2 > 50 ] <- NA

tb_filtered$`C1.整體顏色(原果色)`[is.na(tb_filtered$lyco_1)==F] %>% table()
#咖啡    紅  粉紅  淡黃    黃    綠    橘 
#153 66266  4782  2567  2474  2584  4664 
tb_filtered$`C1.整體顏色(原果色)`[is.na(tb_filtered$lyco_1)] %>% table()
#咖啡   紅 粉紅 淡黃   黃   綠   橘 
#3   7501  673  143   10 2821   39
nolyco <- is.na(tb_filtered$lyco_1)
red <- tb_filtered$`C1.整體顏色(原果色)`=="紅"
pink <- tb_filtered$`C1.整體顏色(原果色)`=="粉紅"
yellow <- tb_filtered$`C1.整體顏色(原果色)`=="黃"
lyellow <- tb_filtered$`C1.整體顏色(原果色)`=="淡黃"
green <- c(1:nrow(tb_filtered)) %in% grep("綠", tb_filtered$`C1.整體顏色(原果色)`)
orange <- tb_filtered$`C1.整體顏色(原果色)`=="橘"
purple <- tb_filtered$`C1.整體顏色(原果色)`=="紫"
brown <- tb_filtered$`C1.整體顏色(原果色)`=="咖啡"
is.na(tb_filtered$`C1.整體顏色(原果色)`) %>% sum() #47
withcol <- is.na(tb_filtered$`C1.整體顏色(原果色)`) == F

tb_brown <- tb_filtered[brown,]
table(tb_brown$`C1.整體顏色(原果色)`)

#Fix ineffective lycopene values for "紅" fruits with median
summary(tb_filtered$lyco_1[red & withcol])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.014  18.889  20.938  21.129  23.312  29.999    7501 
tb_filtered$lyco_1[red & withcol & nolyco] <- mean(tb_filtered$lyco_1[red & withcol], na.rm = T)

#Fix ineffective lycopene values for "粉紅" fruits with median
summary(tb_filtered$lyco_1[pink & withcol])
#Min.   1st Qu.  Median    Mean   3rd Qu.    Max.    NA's 
#0.4719 22.1596 24.2212 24.0971 26.3042   29.9998     673 
tb_filtered$lyco_1[pink & withcol & nolyco] <- mean(tb_filtered$lyco_1[pink & withcol], na.rm = T)

#Fix ineffective lycopene values for "黃" fruits with median
summary(tb_filtered$lyco_1[yellow & withcol])
#Min.   1st Qu.  Median    Mean   3rd Qu.    Max.    NA's 
#0.1733  4.6670  5.8883  5.7993  7.0346   16.0920      10
tb_filtered$lyco_1[yellow & withcol & nolyco] <- mean(tb_filtered$lyco_1[yellow & withcol], na.rm = T)

#Fix ineffective lycopene values for "淡黃" fruits with median
summary(tb_filtered$lyco_1[lyellow & withcol])
#Min.   1st Qu.  Median    Mean  3rd Qu.    Max.    NA's 
#0.00414 1.75824 2.52721 2.58102 3.34332 6.72476     143 
tb_filtered$lyco_1[lyellow & withcol & nolyco] <- mean(tb_filtered$lyco_1[lyellow & withcol], na.rm = T)

#Fix ineffective lycopene values for "綠" fruits with median
summary(tb_filtered$lyco_1[green & withcol])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.0002  1.0064  2.6255  4.8923  6.7061 29.8059    2821 
tb_filtered$lyco_1[green & withcol & nolyco] <- mean(tb_filtered$lyco_1[green & withcol], na.rm = T)

#Fix ineffective lycopene values for "橘" fruits with median
summary(tb_filtered$lyco_1[orange & withcol])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.3299  7.4217  9.8123  9.5770 11.6561 29.5310      39 
tb_filtered$lyco_1[orange & withcol & nolyco] <- mean(tb_filtered$lyco_1[orange & withcol], na.rm = T)

#分離果色含"綠"的紀錄
#GandP_tb <- tb_filtered[grep("綠",tb_filtered$`C1.整體顏色(原果色)`),]
#GandP_tb <- GandP_tb[c(is.na(GandP_tb$lyco_1) | is.na(GandP_tb$lyco_2)) == F,]
#writexl::write_xlsx(GandP_tb, path = "green & purple records.xlsx")

# plots of extimated lycopene------------
lyco_tb <- tb_filtered[,12:17]
lyco_tb <- lyco_tb[is.na(lyco_tb$lyco_1) == F,]
unique(lyco_tb$`C1.整體顏色(原果色)`) %>% sort()
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

# ggplot(data = lyco_tb)+
#   geom_point(mapping = aes(x = lyco_1, y = lyco_2, color = `C1.整體顏色(原果色)`))+
#   scale_color_manual(breaks = c("綠",      
#                                 "綠帶紫","綠帶紫黑","綠紫","紫",
#                                 "咖啡", "紅","粉紅","橘紅","黃紅",
#                                 "橘","黃","淡黃"),
#                      values = c("#6fb510",
#                                 "#087d46","#3b224f","#087d46","#a763dc",
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
sum(WVratio >= 1.2, na.rm = T) #12866
sum(c(WVratio < 1.2 & WVratio > 1.1), na.rm = T) #5932
sum(c(WVratio >= 0.9 & WVratio <= 1.1), na.rm = T) #60413
sum(c(WVratio < 0.9 & WVratio > 0.8), na.rm = T) #13212
sum(WVratio <= 0.8, na.rm = T) #2011
summary(WVratio) #Median = 0.9765
med <- 0.9765

sum(WVratio >= med+0.2, na.rm = T) #14017
sum(c(WVratio < med+0.2 & WVratio > med+0.1), na.rm = T) #6818
sum(c(WVratio >= med-0.1 & WVratio <= med+0.1), na.rm = T) #63674
sum(c(WVratio < med-0.1 & WVratio > med-0.2), na.rm = T) #8623
sum(WVratio <= med-0.2, na.rm = T) #1302

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
  xlim(c(0,200))+
  xlab("Fruit weight(g)")


# Calculate fruit weight per plant (fwpp)--------
a <- c() # a for accession
p <- c() # p for plant count
fwpp <- c() # fwpp for fruit weight per plant
pass <- c() #check if every plant has 3 or more fruits
for (i in comp_phe){
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
for (i in comp_phe){
  row <- tb_filtered$T2.品系代碼 == i
  dt <- tb_filtered[row,]
  sd <- apply(dt[,c(8:11,16)], MARGIN = 2, FUN = sd, na.rm=T)
  sl <- c(sl, mean(unique(dt$`P5.糖度(brix)`), na.rm=T))
  mean <- apply(dt[,c(8:11,16)], MARGIN = 2, FUN = mean, na.rm=T)
  cv <- c(cv, unname(round(sd/mean*100, digits = 2)))
  avg <- c(avg, round(mean,digits = 2))
}

#make a table of coefficient of variations
dim(cv) <- c(5,260)
cv <- data.frame(cv)
colnames(cv) <- comp_phe
rownames(cv) <- c("fruit length", "fruit width", "fruit weight",
                  "sugar level", "estimated lycopene")

#view the distribution of coefficient of variations across accessions
plot(1:ncol(cv),sort(cv[3,]), ylim = c(0,100)) # Fruit weight
plot(1:ncol(cv),sort(cv[5,]), ylim = c(0,100)) # Lycopene content
hist(as.numeric(cv[3,]), main = "Coefficient of Variation of Fruit Weight \nwithin Each Accession (2019 cool)",
     ylab = "Number of Accessions", xlab = "Coefficient of Variation (%)",
     xlim = c(0,max(cv[3,])), breaks = 20)

hist(as.numeric(cv[4,]), main = "Coefficient of Variation of Sugar Level \nwithin Each Accession (2019 cool)",
     ylab = "Number of Accessions", xlab = "Coefficient of Variation (%)",
     xlim = c(0,max(cv[4,])), breaks = 20)

#make a table of trait averages
dim(avg) <- c(5, 260)
avg <- data.frame(avg)
colnames(avg) <- comp_phe
rownames(avg) <- c("fruit length", "fruit width", "fruit weight",
                   "sugar level","estimated lycopene")
avg <- data.frame(t(avg))
avg$serial_num <- as.numeric(rownames(avg))
avg <- avg[,c(6,1:5)]

avg$sugar.level <- sl

#compare the mean fruit weight provided by TARI and the mean value calculated from their raw data
avg_2019 <- readxl::read_xlsx(path = "../phenotype data.xlsx",
                              sheet = "2018-2019")
avg_2019 <- avg_2019[,c(1,12:14,16,15)]
avg_2019 <- avg_2019[avg_2019$serial_num %in% comp_phe,]
#mean difference between fruit weight calculated from raw data and the provided fruit weight
c(avg$fruit.weight - avg_2019$`fruit_weight(g)_2019cool`) %>% mean() #0.0555
cor(avg$fruit.weight, avg_2019$`fruit_weight(g)_2019cool`) #0.9999299
#numbers of accessions that have fruit weight calculated from raw data identical with the provided fruit weight
c(avg$fruit.weight == avg_2019$`fruit_weight(g)_2019cool`) %>% sum() #23
c(round(avg$fruit.weight,digits = 1) == avg_2019$`fruit_weight(g)_2019cool`) %>% sum() #237

c(avg$fruit.width == avg_2019$`fruit_width(cm)_2019cool`) %>% sum() #258
c(avg$fruit.length == avg_2019$`fruit_length(cm)_2019cool`) %>% sum() #258
unique(tb_out$T2.品系代碼) %>% length() #2 accessions 

#calculate yields and compare them with the yields provided by TARI
yield <- c(pc_tb$FWPP)*33000/1000 #fruit weight per plant*33000/1000
plot(yield, avg_2019$`yield(kg/ha)_2019cool`)
cor(yield, avg_2019$`yield(kg/ha)_2019cool`) #0.9807696
sum(yield == avg_2019$`yield(kg/ha)_2019cool`) #4
yield_comp <- data.frame(accession = comp_phe,
                         yield_from_RAW = yield,
                         yield_2019cool_TARI = avg_2019$`yield(kg/ha)_2019cool`,
                         diff = yield - avg_2019$`yield(kg/ha)_2019cool`)
yield_comp$diff_percent <- yield_comp$diff/yield_comp$yield_2019cool_TARI*100

#writexl::write_xlsx(yield_comp, path = "產量比較_108 涼季.xlsx")

avg$yield <- yield_comp$yield_from_RAW

avg_fromRAW_2019cool <- avg

save(avg_fromRAW_2019cool, file = "./avg_fromRAW_2019cool.RData")

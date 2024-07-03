fixName <- function(data_in){
  for(i in 1:nrow(data_in)){
    data_in[i, 'Name'] <- paste0('SAM', strsplit(data_in[i, 'Name'], '\\-')[[1]][3])
  }
  return(data_in)
}

library(dplyr)

pheno <- read.csv('phenotypes_w_salinity_with_xy.csv')
pheno$Year[pheno$Year == '2017h'] <- '2017'
pheno$Year[pheno$Year == '2017w'] <- '2017'
pheno$SLA_leaf_area <- as.numeric(pheno$SLA_leaf_area)
pheno$SLA_leaf_wt <- as.numeric(pheno$SLA_leaf_wt)
pheno$height <- as.numeric(pheno$height)
pheno$clean_wt_kg <- as.numeric(pheno$clean_wt_kg)
pheno$oil <- as.numeric(pheno$oil)
pheno$Year <- as.factor(pheno$Year)
pheno$Treatment <- as.factor(pheno$Treatment)
for(i in 1:nrow(pheno)){
  pheno[i,'Rep'] <- ifelse(pheno[i, 'Treatment'] == 'Control', paste0(pheno[i,'Rep'],'c'), paste0(pheno[i,'Rep'],'s'))
  pheno[i,'Blk'] <- ifelse(pheno[i, 'Treatment'] == 'Control', paste0(pheno[i, 'Blk'],'c'), paste0(pheno[i,'Blk'],'s'))
}

pheno <- pheno[grepl('SAM', pheno$Name),]
pheno <- fixName(pheno)
pheno <- pheno[pheno$Name != 'SAM272',]

table(pheno$survivor_perc)
pheno <- pheno[pheno$survivor_perc >= 10,]
pheno$salt_avg <- NA
for(i in 1:nrow(pheno)){
  pheno[i, 'salt_avg'] <- pheno[i, 'salinity_0_6'] + pheno[i, 'salinity_6_24']
}
pheno$salt_avg <- sapply(pheno$salt_avg, function(x) ifelse(x == 0, 0, round(x/2,2)))
pheno$salt_avg <- round(pheno$salt_avg, 2)
pheno <- pheno[!(is.na(pheno$Year)),]
range(pheno$salt_avg[pheno$Year == '2016'])
mean(pheno$salt_avg[pheno$Year == '2016' & pheno$salt_avg <= 4.0])
mean(pheno$salt_avg[pheno$Year == '2016' & pheno$salt_avg > 4.0])

range(pheno$salt_avg[pheno$Year == '2017'])
mean(pheno$salt_avg[pheno$Year == '2017' & pheno$salt_avg <= 4.0])
mean(pheno$salt_avg[pheno$Year == '2017' & pheno$salt_avg > 4.0])
#########################################################################
#########################################################################
############Flowering################
as.data.frame(colnames(pheno))
f <- pheno[,c(1:3,5,7,9:12,41,21)]
f <- na.omit(f)
fm <- lm(flower~salt_avg, data =f)
out_f <- influencePlot(fm)
resid_df <- data.frame('fitted' = fitted(fm), 'resid' = resid(fm), 'col' = NA, 'salt' = f$salt_avg)
resid_df$stud_resid <- round(resid_df$resid/sd(resid_df$resid),3)
resid_df$abs_stud_resid <- abs(resid_df$stud_resid)
f2 <- f[-which(resid_df$abs_stud_resid >=2.5),]
fm2 <- lm(flower~salt_avg, data = f2)
out_f2 <- influencePlot(fm2)
resid_df2 <- data.frame('fitted' = fitted(fm2), 'resid' = resid(fm2), 'col' = NA, 'salt' = f2$salt_avg)
resid_df2$stud_resid <- round(resid_df2$resid/sd(resid_df2$resid),3)
resid_df2$abs_stud_resid <- abs(resid_df2$stud_resid)
f3 <- f2[-which(resid_df2$abs_stud_resid >=2.5),]
influencePlot(lm(flower~salt_avg, data = f3))
#########################################################################
#########################################################################
#### leaf area ####
as.data.frame(colnames(pheno))
la <- pheno[,c(1:3, 5,7,9:12,41,13)]
la <- na.omit(la)
lam <- lm(SLA_leaf_area~salt_avg, data = la)
out_la <- influencePlot(lam)
resid_df <- data.frame('fitted' = fitted(lam), 'resid' = resid(lam), 'col' = NA, 'salt' = la$salt_avg)
resid_df$stud_resid <- round(resid_df$resid/sd(resid_df$resid),3)
resid_df$abs_stud_resid <- abs(resid_df$stud_resid)
la2 <- la[-which(resid_df$abs_stud_resid >=2.5),]
lam2 <- lm(SLA_leaf_area~salt_avg, data = la2)
out_la2 <- influencePlot(lam2)
#########################################################################
#########################################################################
#### leaf wt ####
as.data.frame(colnames(pheno))
lwt <- pheno[,c(1:3, 5,7,9:12,41,14)]
lwt <- na.omit(lwt)
lwtm <- lm(SLA_leaf_wt~ salt_avg, data = lwt)
out_lwt <- influencePlot(lam)
resid_df <- data.frame('fitted' = fitted(lwtm), 'resid' = resid(lwtm), 'col' = NA, 'salt' = lwt$salt_avg)
resid_df$stud_resid <- round(resid_df$resid/sd(resid_df$resid),3)
resid_df$abs_stud_resid <- abs(resid_df$stud_resid)
lwt2 <- lwt[-which(resid_df$abs_stud_resid >=2.5),]
lwtm2 <- lm(SLA_leaf_wt~salt_avg, data = lwt2)
out_lwt2 <- influencePlot(lwtm2)
lwt3 <- lwt2[-c(as.numeric(tail(row.names(out_lwt2[order(out_lwt2$CookD), ]), 6))),]
#########################################################################
#########################################################################
#### height ####
as.data.frame(colnames(pheno))
ht <- pheno[,c(1:3, 5,7,9:12,41,24)]
ht <- na.omit(ht)
htm <- lm(height~ salt_avg, data = ht)
out_ht <- influencePlot(htm)
resid_df <- data.frame('fitted' = fitted(htm), 'resid' = resid(htm), 'col' = NA, 'salt' = ht$salt_avg)
resid_df$stud_resid <- round(resid_df$resid/sd(resid_df$resid),3)
resid_df$abs_stud_resid <- abs(resid_df$stud_resid)
ht2 <- ht[-which(resid_df$abs_stud_resid >=2.5),]
htm2 <- lm(height~salt_avg, data = ht2)
out_ht2 <- influencePlot(htm2)
resid_df2 <- data.frame('fitted' = fitted(htm2), 'resid' = resid(htm2), 'col' = NA, 'salt' = ht2$salt_avg)
resid_df2$stud_resid <- round(resid_df2$resid/sd(resid_df2$resid),3)
resid_df2$abs_stud_resid <- abs(resid_df2$stud_resid)
ht3 <- ht2[-which(resid_df2$abs_stud_resid >=2.5),]
influencePlot(lm(height~salt_avg, data = ht3))
# ht3 <- ht2[-c(as.numeric(tail(row.names(out_ht2[order(out_ht2$CookD), ]), 6))),]
#########################################################################
#########################################################################
#### oil ####
as.data.frame(colnames(pheno))
oil <- pheno[,c(1:3, 5,7,9:12,41,37)]
oil <- na.omit(oil)
oilm <- lm(oil_perc~ salt_avg, data = oil)
out_oil <- influencePlot(oilm)
resid_df <- data.frame('fitted' = fitted(oilm), 'resid' = resid(oilm), 'col' = NA, 'salt' = oil$salt_avg)
resid_df$stud_resid <- round(resid_df$resid/sd(resid_df$resid),3)
resid_df$abs_stud_resid <- abs(resid_df$stud_resid)
oil2 <- oil[-which(resid_df$abs_stud_resid >=2.5),]
oilm2 <- lm(oil_perc~salt_avg, data = oil2)
out_oil2 <- influencePlot(oilm2)
# resid_df2 <- data.frame('fitted' = fitted(oilm2), 'resid' = resid(oilm2), 'col' = NA, 'salt' = oil2$salt_avg)
# resid_df2$stud_resid <- round(resid_df2$resid/sd(resid_df2$resid),3)
# resid_df2$abs_stud_resid <- abs(resid_df2$stud_resid)
# oil3 <- oil2[-which(resid_df2$abs_stud_resid >=2.5),]
# influencePlot(lm(oil_perc~salt_avg, data = oil3))
oil3 <- oil2[-c(as.numeric(tail(row.names(out_oil2[order(out_oil2$CookD), ]), 6))),]
#########################################################################
#########################################################################
#### yield ####
as.data.frame(colnames(pheno))
yld <- pheno[,c(1:3, 5,7,9:12,41,36)]
yld <- na.omit(yld)
yldm <- lm(clean_wt_kg~ salt_avg, data = yld)
out_yld <- influencePlot(yldm)
resid_df <- data.frame('fitted' = fitted(yldm), 'resid' = resid(yldm), 'col' = NA, 'salt' = yld$salt_avg)
resid_df$stud_resid <- round(resid_df$resid/sd(resid_df$resid),3)
resid_df$abs_stud_resid <- abs(resid_df$stud_resid)
yld2 <- yld[-which(resid_df$abs_stud_resid >=2.5),]
yldm2 <- lm(clean_wt_kg~salt_avg, data = yld2)
out_yld2 <- influencePlot(yldm2)
resid_df2 <- data.frame('fitted' = fitted(yldm2), 'resid' = resid(yldm2), 'col' = NA, 'salt' = yld2$salt_avg)
resid_df2$stud_resid <- round(resid_df2$resid/sd(resid_df2$resid),3)
resid_df2$abs_stud_resid <- abs(resid_df2$stud_resid)
yld3 <- yld2[-which(resid_df2$abs_stud_resid >=2.5),]
yldm3 <- lm(clean_wt_kg~salt_avg, data = yld3)
out_yld3 <- influencePlot(yldm3)
yld4 <- yld3[-c(as.numeric(tail(row.names(out_yld3[order(out_yld3$CookD), ]), 6))),]


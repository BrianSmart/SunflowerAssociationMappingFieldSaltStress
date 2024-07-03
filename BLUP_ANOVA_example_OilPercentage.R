oil_pheno <- read.csv("oil_noOutlier.csv")

#summary info & slope (plasticity)
for(i in unique(oil_pheno$Name)){
  temp_df <- oil_pheno[oil_pheno$Name == i,]
  temp_m <- lm(oil_perc ~ salt_avg, data = temp_df)
  temp_out <- data.frame('Name' = i, 'total_obs' = nrow(temp_df), 'control_obs' = nrow(temp_df[temp_df$Treatment == 'Control',]), 'saline_obs' = nrow(temp_df[temp_df$Treatment == 'Saline',]),
                                                                                      'max_trait' = max(temp_df[,ncol(temp_df)]), 'min_trait' = min(temp_df[,ncol(temp_df)]),
                        'mean_salinity' = round(mean(temp_df$salinity_6_24),4),'max_salinity' = max(temp_df$salinity_6_24), 'min_salinity' = min(temp_df$salinity_6_24),
                        'sd_resid' = round(sd(temp_m$residuals),4),'slope' = round(coef(temp_m)[2],4))
  
  if(i == unique(oil_pheno$Name)[1]){
    pheno2 <- temp_out
  } else {
    pheno2 <- rbind.data.frame(pheno2, temp_out)
  }
}
row.names(pheno2) <- seq(1, nrow(pheno2),1)
#Filter pheno data, 6 or more observations (control (low salinity) and saline (high salinity))
pheno2 <- pheno2[pheno2$total_obs >= 6,]
pheno2 <- pheno2[pheno2$saline_obs > 1,]
#Sorting/ordering 
pheno2$n2 <- gsub('SAM', '', pheno2$Name)
pheno2$n2 <- as.numeric(pheno2$n2)
pheno2 <- pheno2[order(pheno2$n2),]


####################################################################################
####################################################################################
####################################################################################

# Fitting model in asreml 

library(asreml)
library(ASRgenomics)
library(ASRgwas)
oil_pheno <- read.csv("C:/Users/jmcne/OneDrive/Desktop/Work/Salinity/input/oil_noOutlier.csv")

# str(oil_pheno)
oil_pheno$Year <- as.factor(oil_pheno$Year)
oil_pheno$Blk <- as.factor(oil_pheno$Blk)
oil_pheno$Rep <- as.factor(oil_pheno$Rep)
oil_pheno$Name <- as.factor(oil_pheno$Name)
oil_pheno$Treatment <- as.factor(oil_pheno$Treatment)

oil_pheno <- oil_pheno[oil_pheno$Name %in% pheno2$Name,] #pheno2 generated above
length(unique(oil_pheno$Name))
oil_pheno$col <- sapply(oil_pheno$Treatment, function(x) ifelse(x == 'Saline', 'Red','Blue'))
plot(oil_pheno$salt_avg, oil_pheno$oil_perc, col=oil_pheno$col)
abline(v=4, col='black')
#sort in low salinity (oil_c) and high salinity (oil_s)
oil_c <- oil_pheno[oil_pheno$salt_avg <= 4,]
oil_s <- oil_pheno[oil_pheno$salt_avg > 4,]

#fit a model for low salinity
model1c <- asreml(fixed = oil_perc ~ 1 + Year, 
                 random = ~ Name + Year:Name + Year:Rep + Year:Rep:Blk + Year:salt_avg,
                 maxiter = 1000,
                 data = oil_c)
summary(model1c)$varcomp
plot(fitted(model1c), residuals(model1c))

#format and save BLUP for low salinity
tc <- as.data.frame(model1c$coefficients$random)
tc$ef_name <- row.names(tc)
# t$rnum <- seq(1, nrow(t),1)
t2c <- tc[grep('Name', tc$ef_name),]
t2c <- t2c[!(grepl('Year', t2c$ef_name) == TRUE),]
t2c$Name <- sapply(t2c$ef_name, function(x) strsplit(x, 'Name_')[[1]][2])
colnames(t2c)[1] <- 'control'

#fit a model for high salinity
model1s <- asreml(fixed = oil_perc ~ 1, 
                 random = ~ Name + Year:Name + salt_avg,
                 maxiter = 1000,
                 data = oil_s)
summary(model1s)$varcomp
plot(fitted(model1s), residuals(model1s))

#format and save BLUP for high salinity
ts <- as.data.frame(model1s$coefficients$random)
ts$ef_name <- row.names(ts)
# t$rnum <- seq(1, nrow(t),1)
t2s <- ts[grep('Name', ts$ef_name),]
t2s <- t2s[!(grepl('Year', t2s$ef_name) == TRUE),]
t2s$Name <- sapply(t2s$ef_name, function(x) strsplit(x, 'Name_')[[1]][2])
colnames(t2s)[1] <- 'saline'

#format and save results
oil <- merge(t2c, t2s, by='Name')
oil <- oil[,c(2,4)]
oil <- round(oil,3)
plot(oil$control, oil$saline, cex = 2, pch = 16, ylim = c(min(oil),max(oil)), xlim=c(min(oil),max(oil)))
####################################################################################
####################################################################################
####################################################################################

#### Obtain and save variance components and heritability
calc_denoms <- function(data_in, trait_in){ #function calculates harmonic mean of environments (years) and number of plots 
  env_per_geno <- c()
  rep_per_geno <- c()
  for(i in unique(data_in$Name)){
    temp_df <- data_in[data_in$Name == i,]
    if(sum(is.na(temp_df[,trait_in])) != 0){
      temp_df <- temp_df[-which(is.na(temp_df[,trait_in])),]
    }
    if(nrow(temp_df) != 0){
      #harmonic mean environments
      env_per_geno <- c(env_per_geno, 1/length(unique(temp_df$Year)))
      #harmonic mean of reps
      rep_per_geno <- c(rep_per_geno, 1/nrow(temp_df))
    }
  }
  num_geno <- length(unique(data_in$Name))
  out <- data.frame('denom_env'= num_geno/sum(env_per_geno), 'denom_rep' = num_geno/sum(rep_per_geno))
  # out <- data.frame('denom_rep' = num_geno/sum(rep_per_geno))
  out <- round(out,2)
  return(out)
}

#get harmonic mean of env and number of plots for use in heritability
denoms_control <- calc_denoms(oil_c, 'oil_perc')
denoms_saline <- calc_denoms(oil_s, 'oil_perc')
# Variance components for low salinity
cout <- as.data.frame(summary(model1c)$varcomp)
cout$term <- row.names(cout)
row.names(cout) <- seq(1, nrow(cout),1)
cout$term[cout$term == 'units!R'] <- 'error'
cout$term[cout$term == 'Name'] <- 'Genotype'
cout$term[cout$term == 'Year:Name'] <- 'GxE'
cout <- cout[,c(6,1,2,3,4)]
#Variance components for high salinity
sout <- as.data.frame(summary(model1s)$varcomp)
sout$term <- row.names(sout)
row.names(sout) <- seq(1, nrow(sout),1)
sout$term[sout$term == 'units!R'] <- 'error'
sout$term[sout$term == 'Name'] <- 'Genotype'
sout$term[sout$term == 'Year:Name'] <- 'GxE'
sout <- sout[,c(6,1,2,3,4)]

#heritability for low salinity
h_out_low <- data.frame('trait' = 'oil', 'salinity_level' = 'low','genetic_var' = cout$component[which(cout$term == 'Genotype')], 'gxe' = cout$component[which(cout$term == 'GxE')], 'error' = cout$component[which(cout$term == 'error')])
h_out_low[,c(3:ncol(h_out_low))] <- round(h_out_low[,c(3:ncol(h_out_low))],2)
h_out_low <- cbind.data.frame(h_out_low, denoms_control)
h_out_low$h2 <- round(h_out_low$genetic_var/(h_out_low$genetic_var + (h_out_low$gxe/h_out_low$denom_env) + (h_out_low$error/h_out_low$denom_rep)),2)
h_out_low

#heritability for high salinity
h_out_high <- data.frame('trait' = 'oil', 'salinity_level' = 'high','genetic_var' = sout$component[which(sout$term == 'Genotype')], 'gxe' = sout$component[which(sout$term == 'GxE')], 'error' = sout$component[which(sout$term == 'error')])
h_out_high[,c(3:ncol(h_out_high))] <- round(h_out_high[,c(3:ncol(h_out_high))],2)
h_out_high <- cbind.data.frame(h_out_high, denoms_saline)
h_out_high$h2 <- round(h_out_high$genetic_var/(h_out_high$genetic_var + (h_out_high$gxe/h_out_high$denom_env) + (h_out_high$error/h_out_high$denom_rep)),2)
h_out_high

#format and save
h_out <- rbind.data.frame(h_out_low, h_out_high)

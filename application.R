library(ggplot2)
library(ggpubr)
library(dplyr)
library(refund)
library(lme4)
source("~/Documents/CityU/research/Project on Shanghai Actigraph Data/reference/codes/FUI/lfosr3s.R") 
source("~/Documents/CityU/research/Project on Shanghai Actigraph Data/codes/method.R")
## required packages: lme4, refund, dplyr, mgcv, progress, mvtnorm, parallel

setwd("~/Documents/CityU/research/Project on Shanghai Actigraph Data")

sh <- get(load("shanghai.RData"))
sh_all <- get(load("sh_all.RData"))
sh_out <- sh_out[-which(duplicated(sh_out$id)), ]
sh_ori <- read.csv("sh_ori.csv")
colnames(sh_ori)[1] <- "id"
sh_ori <- sh_ori[-which(duplicated(sh_ori$id)), ]

sh_district <- readxl::read_xlsx("上海学校基本信息2.xlsx")  

sh_out <- merge(sh_out, sh_ori[, c(which(colnames(sh_ori) == "学生姓名"|colnames(sh_ori) =="num_class"|colnames(sh_ori) =="school_num"), 
                                   grep("basic", colnames(sh_ori)))],
                by=c("学生姓名", "num_class", "school_num"))

sh_out$bmi <- sh_out$weight/(sh_out$height/100)^2
sh_out$bmi_health <- ifelse(sh_out$bmi <= 18.5, "underweight", 
                            ifelse(sh_out$bmi > 18.5 & sh_out$bmi <= 24.9, "normal", 
                                   ifelse(sh_out$bmi > 24.9 & sh_out$bmi <= 29.9, "overweight", "obese")))
sh_out$bmi_health <- factor(sh_out$bmi_health, levels = c("underweight", "normal", "overweight", "obese"))


#dass vs district
sh_out$district <- ifelse(sh_out$school %in% sh_district$学校[sh_district$区县 == "崇明"], "Chongming",
                          ifelse(sh_out$school %in% sh_district$学校[sh_district$区县 == "嘉定"], "Jiading",
                                 ifelse(sh_out$school %in% sh_district$学校[sh_district$区县 == "金山"], "Jinshan",
                                        ifelse(sh_out$school %in% sh_district$学校[sh_district$区县 == "静安"], "Jingan", 
                                               ifelse(sh_out$school %in% sh_district$学校[sh_district$区县 == "浦东"], "Pudong", "Changning")))))


sh_out$grade_num <- ifelse(sh_out$grade == "七" | sh_out$grade ==  "七年级" | sh_out$grade ==  "初一", "7",
                           ifelse(sh_out$grade == "八" | sh_out$grade ==  "八年级" | sh_out$grade ==  "初二", "8", 
                                  ifelse(sh_out$grade == "九" | sh_out$grade ==  "九年级" | sh_out$grade ==  "初三", "9", 
                                         ifelse(sh_out$grade == "高一", "10", 
                                                ifelse(sh_out$grade == "高二" | sh_out$grade == "高二年级", "11", 
                                                       ifelse(sh_out$grade == "高三", "12", "6"))))))
sh_out$grade_num <- factor(sh_out$grade_num, levels = c(6:12))


#####################
#####################

jx <- get(load("jiangxi.RData"))
jx_all <- get(load("jx_all2.RData"))

jx_out$age <- as.numeric((as.Date(jx_out$startdate)-as.Date(jx_out$birthday))/365.25)

jx_out$district_new <- ifelse(jx_out$district == "信州区", "Xinzhou", 
                              ifelse(jx_out$district == "婺源县", "Wuyuan", 
                                     ifelse(jx_out$district == "玉山县", "Yushan", "Boyang")))
jx_out$grade_num <- ifelse(jx_out$new_grade ==  "初一", "7",
                       ifelse(jx_out$new_grade == "初二", "8", 
                              ifelse(jx_out$new_grade == "初三", "9", 
                                     ifelse(jx_out$new_grade == "高一", "10", 
                                            ifelse(jx_out$new_grade == "高二", "11", "12")))))
jx_out$grade_num <- factor(jx_out$grade_num, levels = 7:12)

jx_out$urban_new <- ifelse(jx_out$urban ==  "农业" | jx_out$urban ==  "农村", "rural",
                           ifelse(jx_out$urban  == "城市" |jx_out$urban  == "城镇" , "urban", NA))
jx_out$judu_new <- ifelse(jx_out$judu ==  "住宿" | jx_out$judu  == "寄宿" | jx_out$judu  == "寄宿生" | 
                            jx_out$judu  == "寄读" | jx_out$judu==  "住校", "lodge",
                                 ifelse(jx_out$judu  == "走读" | jx_out$judu  == "走读生", "extern", NA))
#jx_out <- jx_out[jx_out$judu_new != "unknown",]


SF12v2_gh <- ifelse(jx_out$SF1 == "1", 5,
                           ifelse(jx_out$SF1 == "2", 4.4, 
                                  ifelse(jx_out$SF1 == "3", 3.4, 
                                         ifelse(jx_out$SF1 == "4", 2, 1))))
SF12v2_pf1 <- as.numeric(jx_out$SF2a)
SF12v2_pf2 <- as.numeric(jx_out$SF2b)
SF12v2_rp1 <- as.numeric(jx_out$SF3a)
SF12v2_rp2 <- as.numeric(jx_out$SF3b)
SF12v2_re1 <- as.numeric(jx_out$SF4a)
SF12v2_re2 <- as.numeric(jx_out$SF4b)
SF12v2_bp <- ifelse(jx_out$SF5 == "1", 5,
                           ifelse(jx_out$SF1 == "2", 4, 
                                  ifelse(jx_out$SF1 == "3", 3, 
                                         ifelse(jx_out$SF1 == "4", 2, 1))))
SF12v2_mh1 <- ifelse(jx_out$SF6a == "1", 5,
                          ifelse(jx_out$SF6a == "2", 4, 
                                 ifelse(jx_out$SF6a == "3", 3, 
                                        ifelse(jx_out$SF6a == "4", 2, 1))))
SF12v2_vt <- ifelse(jx_out$SF6b == "1", 5,
                           ifelse(jx_out$SF6b == "2", 4, 
                                  ifelse(jx_out$SF6b == "3", 3, 
                                         ifelse(jx_out$SF6b == "4", 2, 1))))
SF12v2_mh2 <- as.numeric(jx_out$SF6c)
SF12v2_sf <- as.numeric(jx_out$SF7)
SF12v2_pf <- SF12v2_pf1 + SF12v2_pf2
SF12v2_PFtransformed <- SF12v2_tpf <- (SF12v2_pf-2)/4*100

SF12v2_RPraw <- SF12v2_rp <- SF12v2_rp1 + SF12v2_rp2
SF12v2_RPtransformed <- SF12v2_trp <- (SF12v2_rp-2)/8*100
SF12v2_BPtransformed <- SF12v2_tbp <- (SF12v2_bp-1)/4*100
SF12v2_GHtransformed <- SF12v2_tgh <- (SF12v2_gh-1)/4*100
SF12v2_VTtransformed <- SF12v2_tvt <- (SF12v2_vt-1)/4*100
SF12v2_SFtransformed <- SF12v2_tsf <- (SF12v2_sf-1)/4*100

SF12v2_REraw <- SF12v2_re <- SF12v2_re1+SF12v2_re2
SF12v2_REtransformed <- SF12v2_tre <- (SF12v2_re-2)/8*100
SF12v2_MHraw <- SF12v2_mh <- SF12v2_mh1 + SF12v2_mh2
SF12v2_MHtransformed <-  SF12v2_tmh <- (SF12v2_mh-2)/8*100

SF12v2_pf_z<-(SF12v2_tpf-87.6)/22.5
SF12v2_rp_z<-(SF12v2_trp-80.0)/22.9
SF12v2_bp_z<-(SF12v2_tbp-78.1)/25.0
SF12v2_gh_z<-(SF12v2_tgh-48.3)/27.9
SF12v2_vt_z<-(SF12v2_tvt-62.6)/25.5
SF12v2_sf_z<-(SF12v2_tsf-82.0)/23.9
SF12v2_re_z<-(SF12v2_tre-77.4)/21.4
SF12v2_mh_z<-(SF12v2_tmh-69.1)/18.8
SF12v2_mcs_r<-(SF12v2_pf_z*-0.230)+(SF12v2_rp_z*-0.123)+(SF12v2_bp_z*-0.097)+(SF12v2_gh_z*-0.016)+(SF12v2_vt_z*0.235)+(SF12v2_sf_z*0.269)+(SF12v2_re_z*0.434)+(SF12v2_mh_z*0.486)
SF12v2_pcs_r<-(SF12v2_pf_z*0.424)+(SF12v2_rp_z*0.351)+(SF12v2_bp_z*0.318)+(SF12v2_gh_z*0.250)+(SF12v2_vt_z*0.029)+(SF12v2_sf_z*-0.008)+(SF12v2_re_z*-0.192)+(SF12v2_mh_z*-0.221)
SF12v2_mcs<-(SF12v2_mcs_r*10)+50
SF12v2_pcs<-(SF12v2_pcs_r*10)+50
jx_out$SF12v2_mcs <- SF12v2_mcs
jx_out$SF12v2_pcs <- SF12v2_pcs


## Data manipulation
city <- sh
city_out <- sh_out

## smooth accelerometer dat
#city.sm <- lapply(city, function(i) lowess(i, f=.03)$y)
#city <- lapply(city, function(i)  stats::filter(i, sides=2, filter = (rep(1, 60)/60)))

acc <- as.data.frame(do.call("rbind", city))
acc$watch_id <- names(city)

if (sum(city_out$province  == "江西省") > 0) {
  city_out$watch_id <- city_out$jx_startdate
} else {
  city_out$watch_id <- paste(city_out$watch_sn, city_out$startdate,sep="_")
}
df <- merge(city_out, acc, by = "watch_id") #merge two sets based on watch_sn + start_day


ID <- factor(rep(1:nrow(df), each=7))
visit <- rep(1:7, nrow(df))
gender <- factor(rep(df$new_gender, each=7)) #1 male; 2 female
grade <- factor(rep(df$grade_num, each=7)) #from "pre" to "12"
age <- rep(df$age, each=7)
district <- rep(df$district, each=7)
school <- rep(df$school_num, each=7)
class <- rep(df$num_class, each=7)
anxiety <- rep(df$dass_anxiety, each=7)
depression <- rep(df$dass_depression, each=7)
stress <- rep(df$dass_stress, each=7)
bmi <- rep(df$bmi, each=7)
bmi_health <- rep(df$bmi_health, each=7)
bmi_health <- relevel(bmi_health, ref = "normal")
nvalidwedays <- rep(df$nvalidwedays, each=7)
nvalidwkdays <- rep(df$nvalidwkdays, each=7)
days_week <- unlist(lapply(df$startdate, function(x) weekdays(x=as.Date(seq(7), origin=x)-1)))
days_week <- factor(days_week, levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"))
class_by_day <- paste(class, days_week)
class_by_day <- factor(class_by_day, levels = unique(class_by_day))
work <- factor(ifelse(days_week == "Saturday" | days_week == "Sunday", 0, 1))
anx_cat <- rep(df$anx_cat, each=7)
anx_cat <- factor(anx_cat)
#anx_cat <- factor(ifelse(anx_cat != 4, 0, 1))
dep_cat <- rep(df$dep_cat, each=7)
dep_cat <- factor(dep_cat)
dep_cat <- factor(ifelse(dep_cat != 0, 1, 0))
str_cat <-rep(df$str_cat, each=7)
str_cat <- factor(str_cat)
#str_cat <- factor(ifelse(str_cat != 4, 0, 1))
acips <- rep(df$acips, each=7)
erq_cr <- rep(df$erq_cr, each=7)
erq_es <- rep(df$erq_es, each=7)
sf12v2_mcs <- rep(df$SF12v2_mcs, each=7)
sf12v2_pcs <- rep(df$SF12v2_pcs, each=7)

if (sum(city_out$province  == "江西省") > 0) {
  income <- rep(as.numeric(df$basic4), each=7)
  income[which(income >= 9)] <- NA
  income[which(income >= 8)] <- 7
  
  edu <- rep(as.numeric(df$basic2), each=7)
  edu[which(edu >= 11)] <- NA
  edu[which(edu <= 2)] <- 2
  edu[which(edu >= 9)] <- 8
  
  mother.leave <- rep(as.numeric(df$basic7), each=7)
  mother.leave[is.na(mother.leave)] <- NA
} else {
  income <- rep(df$basic4, each=7)
  income[which(income >= 9)] <- NA# income[which(income == 1)] <- 2
  income[which(income >= 8)] <- 7
  
  edu <- rep(df$basic2, each=7)
  edu[which(edu >= 11)] <- NA
  edu[which(edu <= 2)] <- 2
  edu[which(edu >= 9)] <- 9
}


#only in sh
bri <- rep(df$bri, each=7)
tot <- rep(df$tot, each=7)
inh_test <- factor(rep(df$inh_test, each = 7))
shi_test <- factor(rep(df$shi_test, each = 7))
ini_test <- factor(rep(df$ini_test, each = 7))
wor_test <- factor(rep(df$wor_test, each = 7))
pla_test <- factor(rep(df$pla_test, each = 7))
org_test <- factor(rep(df$org_test, each = 7))
mon_test <- factor(rep(df$mon_test, each = 7))
bri_test <- factor(rep(df$bri_test, each = 7))
mi_test <- factor(rep(df$mi_test, each = 7))
tot_test <- factor(rep(df$tot_test, each = 7))

#only in jx 
judu <- factor(rep(df$judu_new, each = 7))
urban <- factor(rep(df$urban_new, each = 7))

vs <- which(colnames(df) == "V1"); ve <- which(colnames(df) == "V10080")
Y <- lapply(as.data.frame(t(df[,vs:ve])), function(x) matrix(x, nrow = 7, ncol = 1440, byrow = TRUE))
Y <- do.call("rbind", Y)

if (sum(city_out$province  == "江西省") > 0) {
  all <- data.frame(ID=ID, visit=visit, gender=gender,
                    grade=as.numeric(grade), days_week=days_week,judu=judu, urban=urban,
                    income=income, edu=edu, mother.leave = mother.leave,
                    anxiety=anxiety, depression=depression, stress=stress,
                    anx_cat=anx_cat, dep_cat=dep_cat, str_cat=str_cat,
                    acips=acips, erq_cr=erq_cr, erq_es,
                    sf12v2_mcs=sf12v2_mcs, sf12v2_pcs=sf12v2_pcs,
                    Y=I(Y)) #construct new dataframe 
  all<- with(all, all[order(ID, days_week),]) #reorder the week day from Monday to Sunday
  all <- na.omit(all)
} else {
  all <- data.frame(ID=ID, visit=visit, gender=gender,
                    nvalidwedays = nvalidwedays, nvalidwkdays = nvalidwkdays,
                    grade=as.numeric(grade),age=age, days_week=days_week,work=work,
                    bmi=bmi, bmi_health=bmi_health,
                    income=income, edu=edu,
                    district=district, school=school, class=class, class_by_day=class_by_day,
                    anxiety=anxiety, depression=depression, stress=stress,
                    anx_cat=anx_cat, dep_cat=dep_cat, str_cat=str_cat,
                    acips=acips, erq_cr=erq_cr, erq_es,
                    bri=bri,  sf12v2_mcs=sf12v2_mcs, sf12v2_pcs=sf12v2_pcs,
                    inh_test=inh_test, shi_test=shi_test,
                    ini_test=ini_test, wor_test=wor_test,
                    pla_test=pla_test, org_test=org_test,
                    mon_test=mon_test, bri_test=bri_test,
                    mi_test=mi_test, tot_test=tot_test,
                    Y=I(Y)) #construct new dataframe 
  all<- with(all, all[order(ID, days_week),]) #reorder the week day from Monday to Sunday
  all <- na.omit(all)
}

#total activity count
all$tac <- apply(all$Y, 1, sum)

####################
#nvalidwkdays and nvalidwedays

if (sum(city_out$province  == "江西省") == 0){
  n <- unique(all$ID[which(all$nvalidwedays + all$nvalidwkdays >= 5)])
  all <- all[all$ID %in% n, ]
}

acc.all <- all
acc.all <- acc.all[acc.all$days_week != "Saturday" & acc.all$days_week != "Sunday",]
acc.all <- acc.all[acc.all$days_week == "Saturday" | acc.all$days_week == "Sunday",]


####total activity count 
pv <- NULL
coef <- NULL
covariate <-  c("acips", "erq_cr", "erq_es", "sf12v2_mcs", "sf12v2_pcs", "dep_cat")
for (i in covariate){
  fml <- as.formula(paste("tac ~ grade + gender + factor(income) + factor(edu) + judu + urban + ", i, "+ (1 | ID)"))
  tac.result <- lmer(formula = fml,
                     data = acc.all, control = lmerControl(optimizer = "bobyqa"))
  tac.sum <- summary(tac.result)
  coef <- rbind(coef, tac.sum$coefficients[,"Estimate"])
  pv <- rbind(pv, pchisq(tac.sum$coefficients[,"t value"]^2, df = 1, lower.tail = F))
}
rownames(coef) <- covariate
rownames(pv) <- covariate
round(coef/10000,2)
round(pv,2)
round(coef[,"acips"]/10000, 2)
round(pv[,"acips"],2)


lm <- lm(formula = depression ~ edu +income , data = acc.all)
summary(lm)


####
I <- 1200
L <- 500
set.seed(123)
st <-sample(unique(all$ID), size = I)
acc.all <- all[all$ID %in% st,]
#acc.all <- all
#acc.all$Y <- acc.all$Y[,181:360]
#acc.all$Y <- acc.all$Y[,361:720]
#acc.all$Y <- acc.all$Y[,721:1080]
#acc.all$Y <- acc.all$Y[,1080:1440]
acc.all$Y <- acc.all$Y[,floor(seq(1, ncol(acc.all$Y), length.out = L))]
save(acc.all, file="example.RData")

acc.all <- acc.all[acc.all$days_week != "Saturday" & acc.all$days_week != "Sunday",]
acc.all <- acc.all[acc.all$days_week == "Saturday" | acc.all$days_week == "Sunday",]

## FUI
fig <- list()
for (i in 1:6){
  covariate <- c("acips", "erq_cr", "erq_es", "sf12v2_mcs", "sf12v2_pcs", "depression")
  fml <- as.formula(paste("Y ~ grade + gender + income + edu +", covariate[i], "+ (1 | ID) + (1 | class_by_day)"))
  raw <- FALSE
  ptm <- proc.time()
  fit_lfosr3s <- our.lfosr3s(formula =fml, data = acc.all, nknots = 30,
                             family = "gaussian", var = TRUE, analytic = TRUE, raw = raw)
  fuitime <- (proc.time() - ptm)[3]
  pic.FUI <- list()
  for (r in 1:6){
    beta.hat.plt <- data.frame(t = seq(1, ncol(fit_lfosr3s$betaHat), length.out = ncol(fit_lfosr3s$betaHat)),
                               beta = fit_lfosr3s$betaHat[r,],
                               lower = fit_lfosr3s$betaHat[r,] - 2*sqrt(diag(fit_lfosr3s$betaHat.var[,,r])),
                               upper = fit_lfosr3s$betaHat[r,] + 2*sqrt(diag(fit_lfosr3s$betaHat.var[,,r])),
                               lower.joint = fit_lfosr3s$betaHat[r,] - fit_lfosr3s$qn[r]*sqrt(diag(fit_lfosr3s$betaHat.var[,,r])),
                               upper.joint = fit_lfosr3s$betaHat[r,] + fit_lfosr3s$qn[r]*sqrt(diag(fit_lfosr3s$betaHat.var[,,r])))
    tit <- c("intercept","grade", "sex","income", "education", covariate[i])
    p.CI <- ggplot() +
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
            geom_ribbon(aes(x = t, ymax = upper.joint, ymin = lower.joint), data = beta.hat.plt, fill = "gray30", alpha = 0.2) +
            geom_ribbon(aes(x = t, ymax = upper, ymin = lower), data = beta.hat.plt, fill = "gray10", alpha = 0.4) +
            geom_line(aes(x = t, y = beta, color = "Estimate"), data = beta.hat.plt, alpha = 1, lty = 5) +
            scale_colour_manual(name="", values=c("Estimate"="blue3")) +
            scale_x_continuous(breaks = seq(1, ncol(fit_lfosr3s$betaHat), length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) +
            labs(title = tit[r], x="Time of Day") + 
            theme(axis.text.x = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)))
    pic.FUI[[r]] <- p.CI
  }
  fig[[i]] <- list(pic.FUI[[2]],pic.FUI[[3]],pic.FUI[[4]],pic.FUI[[5]],pic.FUI[[6]])
}

ggarrange(fig[[6]][[1]], fig[[6]][[2]], fig[[6]][[3]], fig[[6]][[4]],
          fig[[1]][[5]], fig[[2]][[5]], fig[[3]][[5]], fig[[4]][[5]], fig[[5]][[5]], fig[[6]][[5]], 
          ncol=4, nrow = 2, legend = F)
ggarrange(fig[[6]][[1]], fig[[6]][[2]], fig[[6]][[3]], fig[[6]][[4]],
          fig[[1]][[5]], fig[[2]][[5]], fig[[3]][[5]], fig[[4]][[5]], ncol=4, nrow = 2, legend = F)



######################
raw <- FALSE
ptm <- proc.time()
fit_lfosr3s <- our.lfosr3s(formula = Y ~ grade + gender + depression + (1 | ID), data = acc.all, 
                           family = "gaussian",
                           bs = "cr", nknots = 10,
                           var = TRUE, analytic = TRUE, raw = raw)
fuitime <- (proc.time() - ptm)[3]
pic.FUI <- list()
for (r in 1:6){
  beta.hat.plt <- data.frame(t = seq(1, ncol(fit_lfosr3s$betaHat), length.out = ncol(fit_lfosr3s$betaHat)),
                             beta = fit_lfosr3s$betaHat[r,],
                             lower = fit_lfosr3s$betaHat[r,] - 2*sqrt(diag(fit_lfosr3s$betaHat.var[,,r])),
                             upper = fit_lfosr3s$betaHat[r,] + 2*sqrt(diag(fit_lfosr3s$betaHat.var[,,r])),
                             lower.joint = fit_lfosr3s$betaHat[r,] - fit_lfosr3s$qn[r]*sqrt(diag(fit_lfosr3s$betaHat.var[,,r])),
                             upper.joint = fit_lfosr3s$betaHat[r,] + fit_lfosr3s$qn[r]*sqrt(diag(fit_lfosr3s$betaHat.var[,,r])))
  tit <- c("intercept","grade", "sex", "extern", "lodge", "depression")
  p.CI <- ggplot() +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_ribbon(aes(x = t, ymax = upper.joint, ymin = lower.joint), data = beta.hat.plt, fill = "gray30", alpha = 0.2) +
    geom_ribbon(aes(x = t, ymax = upper, ymin = lower), data = beta.hat.plt, fill = "gray10", alpha = 0.4) +
    geom_line(aes(x = t, y = beta, color = "Estimate"), data = beta.hat.plt, alpha = 1, lty = 5) +
    scale_colour_manual(name="", values=c("Estimate"="blue3")) +
    scale_x_continuous(breaks = seq(1, ncol(fit_lfosr3s$betaHat), length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) +
    labs(title = tit[r], x="Time of Day") + 
    theme(axis.text.x = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)))
  pic.FUI[[r]] <- p.CI
}
ggarrange(pic.FUI[[2]], pic.FUI[[3]], pic.FUI[[4]], pic.FUI[[5]], pic.FUI[[6]], 
          ncol=4, nrow = 2, common.legend = T)



#approximate multiple confidence intervals

pic.FUI1 <- list()
for (i in 1:6){
  r <- i #rth predictor
  beta.hat.plt <- data.frame(t = seq(1, ncol(fit_lfosr3s$betaHat), length.out = ncol(fit_lfosr3s$betaHat)),
                              beta = fit_lfosr3s$betaHat[9,] + fit_lfosr3s$betaHat[r+15,],
                              lower = fit_lfosr3s$betaHat[9,] + fit_lfosr3s$betaHat[r+15,] - 2*sqrt(diag(fit_lfosr3s$betaHat.var[,,9])+diag(fit_lfosr3s$betaHat.var[,,r+15])),
                              upper = fit_lfosr3s$betaHat[9,] + fit_lfosr3s$betaHat[r+15,] + 2*sqrt(diag(fit_lfosr3s$betaHat.var[,,9])+diag(fit_lfosr3s$betaHat.var[,,r+15])))
  tit <- c("stress_Tuesday", "stress_Wednesday", "stress_Thursday","stress_Friday","stress_Saturday", "stress_Sunday")
  p.CI <- ggplot() +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_ribbon(aes(x = t, ymax = upper, ymin = lower), data = beta.hat.plt, fill = "gray10", alpha = 0.4) +
    geom_line(aes(x = t, y = beta, color = "Estimate"), data = beta.hat.plt, alpha = 1, lty = 5) +
    scale_colour_manual(name="", values=c("Estimate"="blue3")) +
    scale_x_continuous(breaks = seq(1, ncol(fit_lfosr3s$betaHat), length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) +
    labs(title = tit[i], x="Time of Day") + 
    theme(axis.text.x = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)))
  pic.FUI1[[i]] <- p.CI
}
ggarrange(pic.FUI1[[1]],pic.FUI1[[2]],pic.FUI1[[3]],pic.FUI1[[4]],pic.FUI1[[5]],pic.FUI1[[6]],
          ncol=3, nrow = 2, common.legend = T)




pt.wise <- function(formula, data, family = "gaussian", argvals = NULL, parallel = FALSE, silent = FALSE){
  library(lme4) ## mixed models
  library(refund) ## fpca.face
  library(dplyr) ## organize lapply results
  library(progress) ## display progress bar
  library(mgcv) ## smoothing in step 2
  library(mvtnorm) ## joint CI
  library(parallel) ## mcapply
  
  if(family != "gaussian") analytic <- FALSE ## bootstrap inference for non-Gaussian family
  
  ## Organize the input
  model_formula <- as.character(formula)
  stopifnot(model_formula[1] == "~" & length(model_formula) == 3)
  
  ##########################################################################################
  ## Step 1
  ##########################################################################################
  if(silent == FALSE) print("Step 1: Massive Univariate Mixed Models")
  
  L <- ncol(data[,model_formula[2]]) ## number of observations on the functional domain
  if(is.null(argvals)) argvals <- 1:L
  
  ## function "unimm" fit univariate mixed model at location l
  unimm <- function(l){
    data$Yl <- unclass(data[,model_formula[2]][,l])
    if(family == "gaussian"){
      fit_uni <- suppressMessages(lmer(formula = as.formula(paste0("Yl ~ ", model_formula[3])), 
                                       data = data, control = lmerControl(optimizer = "bobyqa")))
    }else{
      fit_uni <- suppressMessages(glmer(formula = as.formula(paste0("Yl ~ ", model_formula[3])), 
                                        data = data, family = family, control = glmerControl(optimizer = "bobyqa")))
    }
    betaTilde <- lme4::fixef(fit_uni)
    ub <- betaTilde + 1.96*sqrt(diag(vcov(fit_uni)))
    lb <- betaTilde - 1.96*sqrt(diag(vcov(fit_uni)))
    
    #uTilde <- lme4::ranef(fit_uni, drop=TRUE)[[1]]
    yTilde <- predict(fit_uni)
    
    return(list(betaTilde = betaTilde, ub = ub, lb = lb, yTilde = yTilde))
  }
  
  ## fit massive univariate mixed models
  if(parallel == TRUE){
    massmm <- mclapply(argvals, unimm, mc.cores = detectCores() - 1)
  }else{
    massmm <- lapply(argvals, unimm)
  }
  
  ## obtain betaTilde
  betaTilde <- t(lapply(massmm, '[[', 1) %>% bind_rows())
  colnames(betaTilde) <- argvals

  ubTilde <- t(lapply(massmm, '[[', 2) %>% bind_rows())
  lbTilde <- t(lapply(massmm, '[[', 3) %>% bind_rows())
  colnames(ubTilde) <- argvals
  colnames(lbTilde) <- argvals
  
  return(list(betaTilde=betaTilde, ubTilde=ubTilde, lbTilde=lbTilde))
}  



####
acc.all <- all

L <- 200
acc.all$Y <- acc.all$Y[,floor(seq(1, ncol(acc.all$Y), length.out = L))]

acc.all <- acc.all[acc.all$days_week != "Saturday" & acc.all$days_week != "Sunday",]
acc.all <- acc.all[acc.all$days_week == "Saturday" | acc.all$days_week == "Sunday",]

## FUI
fig <- list()
for (i in 1:6){
  print(i)
  covariate <- c("acips", "erq_cr", "erq_es", "sf12v2_mcs", "sf12v2_pcs", "depression")
  fml <- as.formula(paste("Y ~ grade + gender + income + edu +", covariate[i], "+ (1 | ID) + (1 | class_by_day)"))
  ptm <- proc.time()
  fit_lfosr3s <- pt.wise(formula = fml, data = acc.all, 
                         family = "gaussian")
  fuitime <- (proc.time() - ptm)[3]
  pic.FUI <- list()
  for (r in 1:6){
    beta.hat.plt <- data.frame(t = seq(1, ncol(fit_lfosr3s$betaTilde), length.out = ncol(fit_lfosr3s$betaTilde)),
                               beta = fit_lfosr3s$betaTilde[r,],
                               lower = fit_lfosr3s$lbTilde[r,],
                               upper = fit_lfosr3s$ubTilde[r,])
    tit <- c("intercept","grade", "sex", "income", "education", covariate[i])
    p.CI <- ggplot() +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      geom_ribbon(aes(x = t, ymax = upper, ymin = lower), data = beta.hat.plt, fill = "gray10", alpha = 0.4) +
      geom_line(aes(x = t, y = beta, color = "Estimate"), data = beta.hat.plt, alpha = 1, lty = 5, size = 0.8) +
      scale_colour_manual(name="", values=c("Estimate"="blue3")) +
      scale_x_continuous(breaks = seq(1, ncol(fit_lfosr3s$betaTilde), length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) +
      labs(title = tit[r], x="Time of Day") + 
      theme(axis.text.x = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)))
    pic.FUI[[r]] <- p.CI
  }
  fig[[i]] <- list(pic.FUI[[2]],pic.FUI[[3]],pic.FUI[[4]],pic.FUI[[5]],pic.FUI[[6]])
}



ggarrange(fig[[1]][[1]], fig[[1]][[2]], fig[[1]][[3]], fig[[1]][[4]],
          ncol=4, nrow = 2, common.legend = T)


ggarrange(fig[[6]][[1]], fig[[6]][[2]], fig[[6]][[3]], fig[[6]][[4]],
          fig[[1]][[5]] + ylim(-30,40),
          fig[[2]][[5]] + ylim(-40,60), 
          fig[[3]][[5]] + ylim(-60,50), 
          fig[[4]][[5]] + ylim(-35,30), 
          fig[[5]][[5]] + ylim(-30,35),
          fig[[6]][[5]] + ylim(-50,50),
          ncol=4, nrow = 2, legend = F)

ggarrange(fig[[6]][[1]], fig[[6]][[2]], fig[[6]][[3]], fig[[6]][[4]],
          fig[[1]][[5]] + ylim(-15,30),
          fig[[2]][[5]] + ylim(-30,40), 
          fig[[3]][[5]] + ylim(-60,40), 
          fig[[4]][[5]] + ylim(-20,25), 
          fig[[5]][[5]] + ylim(-25,30),
          fig[[6]][[5]] + ylim(-50,50),
          ncol=4, nrow = 2, legend = F)




## FAMM
ptm <- proc.time()
fit_pffr <- pffr(formula = Y ~ age + anx_cat + s(ID, bs = "re"), data = acc.all, family = "gaussian", algorithm = "bam",
                 bs.yindex = list(bs = "ps", k = 15, m = c(2, 1)))
pffrtime <- (proc.time() - ptm)[3]
pffrtime

#For factor variables, a separate smooth function with its own smoothing parameter is estimated for each level of the factor
coef_pffr <- coef(fit_pffr, n1 = L)
betaHat_pffr.var <- betaHat_pffr  <- array(0, dim=c(5, L))

for(i in 1:5){
  if (i == 1) {
    betaHat_pffr[i,] <- as.vector(coef_pffr$smterms[[i]]$value) + coef_pffr$pterms[1]
  } else {
    betaHat_pffr[i,] <- as.vector(coef_pffr$smterms[[i]]$value) 
  }
  betaHat_pffr.var[i,] <- as.vector(coef_pffr$smterms[[i]]$se^2)
}


pic.FAMM <- list() 
for(i in 1:5){
  r <- i #rth predictor
  beta.hat.plt <- data.frame(t = seq(1, L, length.out=L), 
                             beta = betaHat_pffr[r,],
                             upper = betaHat_pffr[r,]+1.96*sqrt(betaHat_pffr.var[r,]),
                             lower = betaHat_pffr[r,]-1.96*sqrt(betaHat_pffr.var[r,]))
  tit <- c("intercept", "age", "anxiety")
  p.CI <- ggplot() +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_ribbon(aes(x = t, ymax = upper, ymin = lower), data = beta.hat.plt, fill = "gray10", alpha = 0.4) +
    geom_line(aes(x = t, y = beta, color = "Estimate"), data = beta.hat.plt, alpha = 1, lty = 5) +
    scale_colour_manual(name="", values=c("Estimate"="blue3")) +
    scale_x_continuous(breaks = seq(1, L,length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) + 
    labs(title = tit[i])
  pic.FAMM[[i]] <- p.CI
}
ggarrange(pic.FAMM[[1]], pic.FAMM[[2]], pic.FAMM[[3]], pic.FAMM[[4]], pic.FAMM[[5]], ncol=3, nrow = 2, common.legend = T)


## Bayesian
ptm <- proc.time()
fit_bayes <- bayes_fosr(formula = Y ~ age + anxiety + re(ID), data = acc.all, Kt = 10, Kp = 10, est.method = "Gibbs", 
                        N.iter = 500, N.burn = 300)
#fit_bayes <- bayes_fosr(formula = Y ~ age + anxiety + re(ID), data = acc.all, Kt = 10, Kp = 10) #variational bayes
bayestime <- (proc.time() - ptm)[3]
bayestime


pic.BAYES<- list() 
for(i in 1:3){
  r <- i #rth predictor
  beta.hat.plt <- data.frame(t = seq(1, L, length.out=L), 
                             beta = fit_bayes$beta.hat[r,],
                             upper = fit_bayes$beta.UB[r,],
                             lower = fit_bayes$beta.LB[r,])
  tit <- c("intercept", "age", "anxiety")
  p.CI <- ggplot() +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_ribbon(aes(x = t, ymax = upper, ymin = lower), data = beta.hat.plt, fill = "gray10", alpha = 0.4) +
    geom_line(aes(x = t, y = beta, color = "Estimate"), data = beta.hat.plt, alpha = 1, lty = 5) +
    scale_colour_manual(name="", values=c("Estimate"="blue3")) +
    scale_x_continuous(breaks = seq(1, L,length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) + 
    labs(title = tit[i])
  pic.BAYES[[i]] <- p.CI
}
ggarrange(pic.BAYES[[1]], pic.BAYES[[2]], pic.BAYES[[3]], ncol=3, nrow = 1, common.legend = T)


## FDboost
library("FDboost")

acc.boost <- list()
acc.boost$Y <- acc.all$Y
acc.boost$age <- acc.all$age
acc.boost$anxiety <- acc.all$anxiety
acc.boost$id <- factor(acc.all$ID)
acc.boost$t <- seq(1:L)


ptm <- proc.time() 
fit_FDboost <- FDboost(formula = Y ~ 1 + bolsc(age, df=3) + bolsc(anxiety, df=3) + brandomc(id, df = 3),
                       data = acc.boost, control = boost_control(mstop = 500), timeformula = ~ bbs(t, df = 4))

bootstrap <- function(data, object, num_bootstrap)
{
  B <- num_bootstrap
  betaHat_boot <- array(NA, dim = c(3, L, B))
  ID.number <- length(unique(data$id))
  for(boots in 1:B){
    sample.ind <- sample(1:ID.number, size = ID.number, replace = TRUE)
    boost.sub <- data
    boost.sub$Y <- boost.sub$Y[boost.sub$id %in% sample.ind, ]
    boost.sub$age <- boost.sub$age[boost.sub$id %in% sample.ind]
    boost.sub$anxiety <- boost.sub$anxiety[boost.sub$id %in% sample.ind]
    boost.sub$id <- boost.sub$id[boost.sub$id %in% sample.ind]
    
    fit_boot <- try(FDboost(formula = as.formula(object$formulaFDboost),
                            data = boost.sub, control = boost_control(mstop = 500), timeformula = as.formula(object$timeformula)))

    coef_boot <- coef(fit_boot, n1 = L, n2 = L)
    betaHat_boot[1,,boots] <- as.vector(coef_boot$smterms[[1]]$value) + coef_boot$offset$value
    betaHat_boot[2,,boots] <- -predict(fit_boot, which=2, newdata=list(age=1, anxiety=1, t=1:L)) 
    betaHat_boot[3,,boots] <- -predict(fit_boot, which=3, newdata=list(age=1, anxiety=1, t=1:L))
  }
  ## obtain bootstrap variance
  betaHat.var <- array(NA, dim = c(L,L,3))
  for(r in 1:3){
    betaHat.var[,,r] <- 1.2*var(t(betaHat_boot[r,,])) ## account for within-subject correlation
  }
  return(betaHat.var)
}
bound <- bootstrap(data = acc.boost, object = fit_FDboost, num_bootstrap = 10)

FDboosttime <- (proc.time() - ptm)[3]
FDboosttime


betaHat_FDboost  <- array(0, dim=c(3, L))

coef_FDboost <- coef(fit_FDboost, n1 = L, n2 = L)
betaHat_FDboost[1,] <- as.vector(coef_FDboost$smterms[[1]]$value) + coef_FDboost$offset$value
betaHat_FDboost[2,] <- -predict(fit_FDboost, which=2, newdata=list(age=1, t=1:L)) 
betaHat_FDboost[3,] <- -predict(fit_FDboost, which=3, newdata=list(anxiety=1, t=1:L))

FDboostbeta.LB <- FDboostbeta.UB <- array(0, dim=c(3, L))
for (r in 1:3){
  FDboostbeta.LB[r,] <- betaHat_FDboost[r,] - 1.96*sqrt(diag(bound[,,r]))
  FDboostbeta.UB[r,] <- betaHat_FDboost[r,] + 1.96*sqrt(diag(bound[,,r]))
}


#betaHat_FDboost[2,] <- apply(coef_FDboost$smterms[[2]]$value, 2, mean)
#betaHat_FDboost[3,] <- apply(coef_FDboost$smterms[[3]]$value, 2, mean)


pic.FDboost <- list() 
for(r in 1:3){
  beta.hat.plt <- data.frame(t = seq(1, L, length.out=L), 
                             beta = betaHat_FDboost[r,],
                             upper = FDboostbeta.UB[r,],
                             lower = FDboostbeta.LB[r,])
  tit <- c("intercept", "age", "anxiety")
  p.CI <- ggplot() +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_ribbon(aes(x = t, ymax = upper, ymin = lower), data = beta.hat.plt, fill = "gray10", alpha = 0.4) +
    geom_line(aes(x = t, y = beta, color = "Estimate"), data = beta.hat.plt, alpha = 1, lty = 5) +
    scale_colour_manual(name="", values=c("Estimate"="blue3")) +
    scale_x_continuous(breaks = seq(1, L,length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) + 
    labs(title = tit[r])
  pic.FDboost[[r]] <- p.CI
}
ggarrange(pic.FDboost[[1]], pic.FDboost[[2]], pic.FDboost[[3]], ncol=3, nrow = 1, common.legend = T)



## plots 
## reorder the signal from Monday to Sunday 
all.reorder <- with(all, all[order(ID, days_week),])
Ynew <- NULL
for (i in 1:length(unique(all.reorder$ID))) {
 new  <- as.vector(t(all.reorder$Y[((i-1)*7+1):(i*7),]))
 Ynew <- rbind(new, Ynew) #new reordered signals
}


plot.subject <- data.frame(y = Ynew[1,], t=rep(1:1440, 7), day = rep(1:7, each=1440)) 

p1 <- ggplot(data = plot.subject, aes(x=1:(7*1440), y=y)) +
  geom_line() +
  labs(x="", y="Accelerometer", title = "Subject 1" ) + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16, margin = unit(c(0, 5, 0, 0), "mm")),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5,margin = unit(c(0, 0, 5, 0), "mm"))) + 
  scale_x_continuous(breaks = seq(720, 9360, length.out = 7), labels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) 


# activity averaged over weekday per age and gender curves
sub.avg <- NULL
for (i in 0:1){
  sub <- all[all$days_week != "Saturday" | all$days_week != "Sunday", ]
  sub <- sub$Y[sub$ID %in% unique(sub$ID[which(sub$dep_cat==i)]),]
  sub.avg <-c(sub.avg, apply(sub, 2, mean))
}
sub.avg <- data.frame(time=rep(1:1440,2), avg = sub.avg, cat = rep(paste0("Depression",0:1), each=1440))

p2 <- ggplot(data = sub.avg, aes(x=time, y=avg,  color=cat)) +
  geom_line(alpha=0.6) +
  labs(x="", y="Accelerometer", title = "Weekday Average" ) + 
  theme_bw() +
  ylim(0,4500) +
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16, margin = unit(c(0, 5, 0, 0), "mm")),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5,margin = unit(c(0, 0, 5, 0), "mm")), 
        legend.position="bottom", legend.title = element_blank(),  legend.text = element_text(size=15) ) + 
  scale_x_continuous(breaks = seq(1, 1440,length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) 

sub.avg <- NULL
for (i in 0:1){
  sub <- all[all$days_week == "Saturday" | all$days_week == "Sunday", ]
  sub <- sub$Y[sub$ID %in% unique(sub$ID[which(sub$dep_cat==i)]),]
  sub.avg <-c(sub.avg, apply(sub, 2, mean))
}
sub.avg <- data.frame(time=rep(1:1440,2), avg = sub.avg, cat = rep(paste0("Depression",0:1), each=1440))


p3 <- ggplot(data = sub.avg, aes(x=time, y=avg,  color=cat)) +
  geom_line(alpha=0.6) +
  labs(x="", y="Accelerometer", title = "Weekend Average" ) + 
  theme_bw() +
  ylim(0,4500) +
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16, margin = unit(c(0, 5, 0, 0), "mm")),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5,margin = unit(c(0, 0, 5, 0), "mm")), 
        legend.position="bottom", legend.title = element_blank(),  legend.text = element_text(size=15) ) + 
  scale_x_continuous(breaks = seq(1, 1440,length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep=""))

ggarrange(p1, 
          ggarrange(p2, p3, common.legend = T, legend = "bottom"), nrow = 2, align = "h")


##TUGoverSex

sub.avg <- NULL
for (i in 1:2){
  sub <- all[all$days_week != "Saturday" & all$days_week != "Sunday", ]
  sub <- sub$Y[sub$ID %in% unique(sub$ID[which(sub$gender==i)]),]
  sub.avg <-c(sub.avg, apply(sub, 2, mean))
}

sub.avg <- data.frame(time=rep(1:1440,2), avg = sub.avg, cat = rep(1:2, each=1440))

p1 <- ggplot(data = sub.avg, aes(x=time, y=avg,  color=factor(cat))) +
  geom_line(alpha=0.9) +
  labs(x="", y="Accelerometer", title = "Weekday Average" ) + 
  theme_bw() +
  ylim(0,4500) +
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16, margin = unit(c(0, 5, 0, 0), "mm")),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5,margin = unit(c(0, 0, 5, 0), "mm")), 
        legend.title = element_text(size=15), legend.text = element_text(size=15) ) + 
  scale_x_continuous(breaks = seq(1, 1440,length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) + 
  scale_colour_discrete(name="Sex", labels=c("Male","Female"))


sub.avg <- NULL
for (i in 1:2){
  sub <- all[all$days_week == "Saturday" | all$days_week == "Sunday", ]
  sub <- sub$Y[sub$ID %in% unique(sub$ID[which(sub$gender==i)]),]
  sub.avg <-c(sub.avg, apply(sub, 2, mean))
}

sub.avg <- data.frame(time=rep(1:1440,2), avg = sub.avg, cat = rep(1:2, each=1440))

p2 <- ggplot(data = sub.avg, aes(x=time, y=avg,  color=factor(cat))) +
  geom_line(alpha=0.9) +
  labs(x="", y="Accelerometer", title = "Weekend Average" ) + 
  theme_bw() +
  ylim(0,4500) +
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16, margin = unit(c(0, 5, 0, 0), "mm")),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5,margin = unit(c(0, 0, 5, 0), "mm")), 
        legend.title = element_text(size=15), legend.text = element_text(size=15) ) + 
  scale_x_continuous(breaks = seq(1, 1440,length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) + 
  scale_colour_discrete(name="Sex", labels=c("Male","Female"))

ggarrange(p1, p2, common.legend = T, legend = "right")



##TUGoverGrade

sub.avg <- NULL
for (i in 1:7){
  sub <- all[all$days_week != "Saturday" & all$days_week != "Sunday", ]
  sub <- sub$Y[sub$ID %in% unique(sub$ID[which(sub$grade==i)]),]
  avg <- apply(sub, 2, mean)
  avg <- stats::filter(avg, sides=2, filter = (rep(1, 20)/20))
  sub.avg <-c(sub.avg, avg)
}

sub.avg <- data.frame(time=rep(1:1440,7), avg = sub.avg, cat = rep(1:7, each=1440))

p1 <- ggplot(data = sub.avg, aes(x=time, y=avg,  color=factor(cat))) +
  geom_line(alpha=0.9) +
  labs(x="", y="Accelerometer", title = "Weekday Average" ) + 
  theme_bw() +
  ylim(0,4500) +
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16, margin = unit(c(0, 5, 0, 0), "mm")),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5,margin = unit(c(0, 0, 5, 0), "mm")), 
        legend.title = element_text(size=15), legend.text = element_text(size=15) ) + 
  scale_x_continuous(breaks = seq(1, 1440,length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) + 
  scale_colour_discrete(name="Grade", labels=6:12)



sub.avg <- NULL
for (i in 1:7){
  sub <- all[all$days_week == "Saturday" | all$days_week == "Sunday", ]
  sub <- sub$Y[sub$ID %in% unique(sub$ID[which(sub$grade==i)]),]
  avg <- apply(sub, 2, mean)
  avg <- stats::filter(avg, sides=2, filter = (rep(1, 20)/20))
  sub.avg <-c(sub.avg, avg)
}

sub.avg <- data.frame(time=rep(1:1440,7), avg = sub.avg, cat = rep(1:7, each=1440))

p2 <-  ggplot(data = sub.avg, aes(x=time, y=avg,  color=factor(cat))) +
  geom_line(alpha=0.9) +
  labs(x="", y="Accelerometer", title = "Weekday Average" ) + 
  theme_bw() +
  ylim(0,4500) +
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16, margin = unit(c(0, 5, 0, 0), "mm")),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5,margin = unit(c(0, 0, 5, 0), "mm")), 
        legend.title = element_text(size=15), legend.text = element_text(size=15) ) + 
  scale_x_continuous(breaks = seq(1, 1440,length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) + 
  scale_colour_discrete(name="Grade", labels=6:12)

ggarrange(p1, p2, common.legend = T, legend = "right")



##TUGoverEdu

sub.avg <- NULL
for (i in 2:9){
  sub <- all[all$days_week != "Saturday" & all$days_week != "Sunday", ]
  sub <- sub$Y[sub$ID %in% unique(sub$ID[which(sub$edu==i)]),]
  avg <- apply(sub, 2, mean)
  avg <- stats::filter(avg, sides=2, filter = (rep(1, 20)/20))
  sub.avg <-c(sub.avg, avg)
}

sub.avg <- data.frame(time=rep(1:1440,8), avg = sub.avg, cat = rep(1:8, each=1440))

p1 <- ggplot(data = sub.avg, aes(x=time, y=avg,  color=factor(cat))) +
  geom_line(alpha=0.9) +
  labs(x="", y="Accelerometer", title = "Weekday Average" ) + 
  theme_bw() +
  ylim(0,4500) +
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16, margin = unit(c(0, 5, 0, 0), "mm")),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5,margin = unit(c(0, 0, 5, 0), "mm")), 
        legend.title = element_text(size=13), legend.text = element_text(size=13), legend.text.align = 0 ) + 
  scale_x_continuous(breaks = seq(1, 1440,length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) + 
  scale_colour_discrete(name="Education", labels=c(expression(""<"Primary school"), "Primary school", "Middle school", "High school", 
                                                   "Technical secondary school", "College", "Bachelor", expression("">="Master")))


sub.avg <- NULL
for (i in 2:9){
  sub <- all[all$days_week == "Saturday" | all$days_week == "Sunday", ]
  sub <- sub$Y[sub$ID %in% unique(sub$ID[which(sub$edu==i)]),]
  avg <- apply(sub, 2, mean)
  avg <- stats::filter(avg, sides=2, filter = (rep(1, 20)/20))
  sub.avg <-c(sub.avg, avg)
}

sub.avg <- data.frame(time=rep(1:1440,8), avg = sub.avg, cat = rep(1:8, each=1440))

p2 <-  ggplot(data = sub.avg, aes(x=time, y=avg,  color=factor(cat))) +
  geom_line(alpha=0.9) +
  labs(x="", y="Accelerometer", title = "Weekday Average" ) + 
  theme_bw() +
  ylim(0,4500) +
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16, margin = unit(c(0, 5, 0, 0), "mm")),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5,margin = unit(c(0, 0, 5, 0), "mm")), 
        legend.title = element_text(size=13), legend.text = element_text(size=13) ) + 
  scale_x_continuous(breaks = seq(1, 1440,length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) + 
  scale_colour_discrete(name="Education", labels=c(expression(""<"Primary school"), "Primary school", "Middle school", "High school", 
                                                   "Technical secondary school", "College", "Bachelor", expression("">="Master")))


ggarrange(p1, p2, common.legend = T, legend = "right")



##TUGoverIncome

sub.avg <- NULL
for (i in 1:7){
  sub <- all[all$days_week != "Saturday" & all$days_week != "Sunday", ]
  sub <- sub$Y[sub$ID %in% unique(sub$ID[which(sub$income==i)]),]
  avg <- apply(sub, 2, mean)
  avg <- stats::filter(avg, sides=2, filter = (rep(1, 20)/20))
  sub.avg <-c(sub.avg, avg)
}

sub.avg <- data.frame(time=rep(1:1440,7), avg = sub.avg, cat = rep(1:7, each=1440))

p1 <- ggplot(data = sub.avg, aes(x=time, y=avg,  color=factor(cat))) +
  geom_line(alpha=0.9) +
  labs(x="", y="Accelerometer", title = "Weekday Average" ) + 
  theme_bw() +
  ylim(0,4500) +
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16, margin = unit(c(0, 5, 0, 0), "mm")),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5,margin = unit(c(0, 0, 5, 0), "mm")), 
        legend.title = element_text(size=13), legend.text = element_text(size=13), legend.text.align = 0 ) + 
  scale_x_continuous(breaks = seq(1, 1440,length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) + 
  scale_colour_discrete(name="Income", labels=c(expression(""<="10k"), "10-30k", "30-50k", "50-100k", 
                                        "100-150k", "150-300k", expression("">="300k"))) 

sub.avg <- NULL
for (i in 1:7){
  sub <- all[all$days_week == "Saturday" | all$days_week == "Sunday", ]
  sub <- sub$Y[sub$ID %in% unique(sub$ID[which(sub$income==i)]),]
  avg <- apply(sub, 2, mean)
  avg <- stats::filter(avg, sides=2, filter = (rep(1, 20)/20))
  sub.avg <-c(sub.avg, avg)
}

sub.avg <- data.frame(time=rep(1:1440,7), avg = sub.avg, cat = rep(1:7, each=1440))


p2 <-  ggplot(data = sub.avg, aes(x=time, y=avg,  color=factor(cat))) +
  geom_line(alpha=0.9) +
  labs(x="", y="Accelerometer", title = "Weekday Average" ) + 
  theme_bw() +
  ylim(0,4500) +
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16, margin = unit(c(0, 5, 0, 0), "mm")),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5,margin = unit(c(0, 0, 5, 0), "mm")), 
        legend.title = element_text(size=13), legend.text = element_text(size=13) ) + 
  scale_x_continuous(breaks = seq(1, 1440,length.out = 5), labels = paste(c(0,6,12,18,24), ":00", sep="")) + 
  scale_colour_discrete(name="Income", labels=c(expression(""<="10k"), "10-30k", "30-50k", "50-100k", 
                                        "100-150k", "150-300k", expression("">="300k")))
ggarrange(p1, p2, common.legend = T, legend = "right")



#####MVPA, LA. SED
Ymvp <- ifelse(Ynew>=3600, "MVPA", ifelse(Ynew<200, "SED", "LA"))

Ymvp.sum <- Ymvp[,as.vector(sapply(1:5, function(x) x*(1081:1440)))] #weekday
Ymvp.sum <- base::apply(Ymvp.sum, 1, function(x) table(factor(x, levels = c("SED", "LA", "MVPA")))/5)
Ymvp.sum <- t(Ymvp.sum)
#Ymvp.sum <- do.call("rbind", Ymvp.sum)
Ymvp.sum <- cbind(Ymvp.sum, grade=sapply(1:nrow(Ynew), function(x) unique(all.reorder$grade[all.reorder$ID==x])))

Ymvp.sum.plot <- data.frame(Ymvp.sum) %>% dplyr::group_by(grade) %>%
  dplyr::summarise(MVPAm=mean(MVPA), LAm=mean(LA), SEDm=mean(SED))

Ymvp.sum.plot <- reshape(Ymvp.sum.plot, idvar="grade",
                         varying=c("MVPAm", "LAm", "SEDm"), 
                         times=c("MVPA", "LA", "SED"), 
                         new.row.names = 1:1000, v.names=c("avg"), direction="long")

p1 <- ggplot(Ymvp.sum.plot, aes(x=grade, y=avg, color=time)) + 
  geom_point(size=1.3) + 
  geom_line(size=1.3) +
  theme_bw() + 
  scale_color_discrete(name="") +
  labs(y = "Minutes per day", x = "Grade", title = "Weekday") + 
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14),
        legend.text = element_text(size=12), legend.title = element_text(size=14), 
        plot.title = element_text(size=20)) + 
  scale_x_continuous(breaks=1:7, labels=6:12) + 
  coord_cartesian(ylim = c(0, 300))


Ymvp.sum <- Ymvp[,as.vector(sapply(6:7, function(x) x*(1081:1440)))] #weekend
Ymvp.sum <- base::apply(Ymvp.sum, 1, function(x) table(factor(x, levels = c("SED", "LA", "MVPA")))/2)
Ymvp.sum <- t(Ymvp.sum)
#Ymvp.sum <- do.call("rbind", Ymvp.sum)
Ymvp.sum <- cbind(Ymvp.sum, grade=sapply(1:nrow(Ynew), function(x) unique(all.reorder$grade[all.reorder$ID==x])))

Ymvp.sum.plot <- data.frame(Ymvp.sum) %>% dplyr::group_by(grade) %>%
  dplyr::summarise(MVPAm=mean(MVPA), LAm=mean(LA), SEDm=mean(SED))

Ymvp.sum.plot <- reshape(Ymvp.sum.plot, idvar="grade",
        varying=c("MVPAm", "LAm", "SEDm"), 
        times=c("MVPA", "LA", "SED"), 
        new.row.names = 1:1000, v.names=c("avg"), direction="long")

p2 <- ggplot(Ymvp.sum.plot, aes(x=grade, y=avg, color=time)) + 
  geom_point(size=1.3) + 
  geom_line(size=1.3) +
  theme_bw() + 
  scale_color_discrete(name="") +
  labs(y = "Minutes per day", x = "Grade", title = "Weekend") + 
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14),
        legend.text = element_text(size=12), legend.title = element_text(size=14),
        plot.title = element_text(size=20)) + 
  scale_x_continuous(breaks=1:7, labels=6:12) + 
  coord_cartesian(ylim = c(0, 300))

ggarrange(p1, p2, common.legend = T, legend = "right")





par(mfrow=c(1,2))
t <- 1:(7*1440)
plot(plot.subject.sm$y~t, type="l", ylab = "Accelerometer", xlab = "", main = "Subject1",  xaxt = "n")
axis(1, at=seq(720, 9360, length.out = 7), labels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"))

# prepare variables.
x <- 1:1440
y <- 1:7
z <- t(all.reorder[all.reorder$ID == 2, "Y"])
z <- array(city.sm$`MOS2D30170805_2017-12-15`, dim = c(1440,7))

# plot the 3D surface
 persp(x, y, z.sm,
      main = "Subject 1",
      xlab = "Time of Day",
      ylab = "Week Day",
      zlab = "Accelerometer",
      theta = 50, phi = 35,
      col = "yellow", shade = 0.005)


heatmap(z.sm, Rowv = NA, Colv = "Rowv")



####################






                                
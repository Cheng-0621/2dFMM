## Data manipulation
city <- jx
city_out <- jx_out

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
  all <- na.omit(all)
} else {
  all <- data.frame(ID=ID, visit=visit, gender=gender,
                    nvalidwedays = nvalidwedays, nvalidwkdays = nvalidwkdays,
                    grade=as.numeric(grade),age=age, days_week=days_week,work=work,
                    bmi=bmi, bmi_health=bmi_health,
                    income=income, edu=edu,
                    district=district, school=school,
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
weekday <- acc.all[acc.all$days_week != "Saturday" & acc.all$days_week != "Sunday",]
weekend <- acc.all[acc.all$days_week == "Saturday" | acc.all$days_week == "Sunday",]


#######Shanghai##########
if (sum(weekday$urban == "urban") > 0){
  data.tac <- weekday %>% group_by(ID) %>% 
    dplyr::summarise(tac.mean = mean(tac, na.rm = T), 
                     grade = unique(grade), 
                     gender = unique(gender),
                     income = unique(income), 
                     edu = unique(edu),
                     acips.mean = mean(acips, na.rm = T),
                     erq_cr.mean = mean(erq_cr, na.rm = T), 
                     erq_es.mean = mean(erq_es, na.rm = T),
                     sf12v2_mcs.mean = mean(sf12v2_mcs, na.rm = T),
                     sf12v2_pcs.mean = mean(sf12v2_pcs, na.rm = T),
                     dep.mean = mean(depression, na.rm = T), 
                     school = unique(school), 
                     judu = unique(judu), 
                     mother.leave = unique(mother.leave), 
                     urban = unique(urban))
} else {
  data.tac <- weekday %>% group_by(ID) %>% 
    dplyr::summarise(tac.mean = mean(tac, na.rm = T), 
                     grade = unique(grade), 
                     gender = unique(gender),
                     income = unique(income), 
                     edu = unique(edu),
                     acips.mean = mean(acips, na.rm = T),
                     erq_cr.mean = mean(erq_cr, na.rm = T), 
                     erq_es.mean = mean(erq_es, na.rm = T),
                     sf12v2_mcs.mean = mean(sf12v2_mcs, na.rm = T),
                     sf12v2_pcs.mean = mean(sf12v2_pcs, na.rm = T),
                     dep.mean = mean(depression, na.rm = T), 
                     school = unique(school))
}

give.n <- function(x){
  return(data.frame(y = quantile(x + 0.3, probs = 0.99), label = length(x)))
                    #label = paste0("n = ",length(x))))
}


p1 <- ggplot(data.tac, aes(x=as.factor(gender), y=log(tac.mean), fill = as.factor(gender))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA, width = 0.5) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 4) +
  scale_fill_discrete(name = "Sex", labels = c("Male", "Female")) + 
  #guides(fill=guide_legend(title="Sex")) + 
  labs(y = "log(Total Activity Count) \n Weekday", title = "Gender") + 
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12),  axis.text.x = element_text(size=12), 
        legend.position="none") + 
  scale_x_discrete(breaks=1:2, labels=c("Male", "Female")) + 
  coord_cartesian(ylim = c(13, 16))

p2 <- ggplot(data.tac, aes(x=factor(grade), y=log(tac.mean), fill = as.factor(grade))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 3.5) +
  scale_fill_discrete(name = "Grade", labels = 7:12) + 
  #guides(fill=guide_legend(title="Grade")) + 
  labs(title = "Grade") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), 
        legend.position="none") + 
  scale_x_discrete(breaks=1:7, labels=6:12) + 
  coord_cartesian(ylim = c(13, 16))

p3 <- ggplot(subset(data.tac, !is.na(income)), aes(x=factor(income),  y=log(tac.mean), fill = as.factor(income))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 3.5) +
  scale_fill_discrete(name = "Income", labels = 1:8) + 
  #guides(fill=guide_legend(title="Income")) + 
  labs(title = "Income (¥)") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=10, angle = 30, hjust=1), 
        legend.position="none") +  
  scale_x_discrete(breaks=1:7, labels=c(expression(""<="10k"), "10-30k", "30-50k", "50-100k", 
                                        "100-150k", "150-300k", expression("">="300k"))) + 
  coord_cartesian(ylim = c(13, 16))

p4 <- ggplot(subset(data.tac, !is.na(edu)), aes(x=factor(edu), y=log(tac.mean), fill = as.factor(edu))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 3.5) +
  scale_fill_discrete(name = "Education", labels = 1:9) + 
  #guides(fill=guide_legend(title="Education")) +
  labs(title = "Education") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=10, angle = 30, hjust=1), 
        legend.position="none") +  
  scale_x_discrete(breaks=2:9, labels=c(expression(""<"Primary school"), "Primary school", "Middle school", "High school", 
                                        "Technical secondary school", "College", "Bachelor", expression("">="Master")))  +
  coord_cartesian(ylim = c(13, 16))
  

ggarrange(p1, p2, p3, p4, 
          ncol=4, nrow=2, align = "h")


if (sum(weekend$urban == "urban") > 0){
  data.tac <- weekend %>% group_by(ID) %>% 
    dplyr::summarise(tac.mean = mean(tac, na.rm = T), 
                     grade = unique(grade), 
                     gender = unique(gender),
                     income = unique(income), 
                     edu = unique(edu),
                     acips.mean = mean(acips, na.rm = T),
                     erq_cr.mean = mean(erq_cr, na.rm = T), 
                     erq_es.mean = mean(erq_es, na.rm = T),
                     sf12v2_mcs.mean = mean(sf12v2_mcs, na.rm = T),
                     sf12v2_pcs.mean = mean(sf12v2_pcs, na.rm = T),
                     dep.mean = mean(depression, na.rm = T), 
                     school = unique(school), 
                     judu = unique(judu), 
                     mother.leave = unique(mother.leave), 
                     urban = unique(urban))
} else {
  data.tac <- weekend %>% group_by(ID) %>% 
    dplyr::summarise(tac.mean = mean(tac, na.rm = T), 
                     grade = unique(grade), 
                     gender = unique(gender),
                     income = unique(income), 
                     edu = unique(edu),
                     acips.mean = mean(acips, na.rm = T),
                     erq_cr.mean = mean(erq_cr, na.rm = T), 
                     erq_es.mean = mean(erq_es, na.rm = T),
                     sf12v2_mcs.mean = mean(sf12v2_mcs, na.rm = T),
                     sf12v2_pcs.mean = mean(sf12v2_pcs, na.rm = T),
                     dep.mean = mean(depression, na.rm = T), 
                     school = unique(school))
}


p8 <- ggplot(data.tac, aes(x=as.factor(gender), y=log(tac.mean), fill = as.factor(gender))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA, width = 0.5) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 4) +
  scale_fill_discrete(name = "Sex", labels = c("Male", "Female")) + 
  #guides(fill=guide_legend(title="Sex")) + 
  labs(y = "log(Total Activity Count) \n Weekday", title = "Gender") + 
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_blank(), 
        axis.text=element_text(size=12), 
        legend.position="none") + 
  scale_x_discrete(breaks=1:2, labels=c("Male", "Female")) + 
  coord_cartesian(ylim = c(13, 16))

p9 <- ggplot(data.tac, aes(x=factor(grade), y=log(tac.mean), fill = as.factor(grade))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 3.5) +
  scale_fill_discrete(name = "Grade", labels = 7:12) + 
  #guides(fill=guide_legend(title="Grade")) + 
  labs(title = "Grade") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text=element_text(size=12), 
        legend.position="none") + 
  scale_x_discrete(breaks=1:7, labels=6:12) + 
  coord_cartesian(ylim = c(13, 16))

p10 <- ggplot(subset(data.tac, !is.na(income)), aes(x=factor(income),  y=log(tac.mean), fill = as.factor(income))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 3.5) +
  scale_fill_discrete(name = "Income", labels = 1:8) + 
  #guides(fill=guide_legend(title="Income")) + 
  labs(title = "Income (¥)") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=10, angle = 30, hjust=1), 
        legend.position="none") +  
  scale_x_discrete(breaks=1:7, labels=c(expression(""<="10k"), "10-30k", "30-50k", "50-100k", 
                                        "100-150k", "150-300k", expression("">="300k"))) + 
  coord_cartesian(ylim = c(13, 16))

p11 <- ggplot(subset(data.tac, !is.na(edu)), aes(x=factor(edu), y=log(tac.mean), fill = as.factor(edu))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 3.5) +
  scale_fill_discrete(name = "Education", labels = 1:9) + 
  #guides(fill=guide_legend(title="Education")) +
  labs(title = "Education") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12),  axis.text.x = element_text(size=10, angle = 30, hjust=1), 
        legend.position="none") +  
  scale_x_discrete(breaks=2:9, labels=c(expression(""<"Primary school"), "Primary school", "Middle school", "High school", 
                                        "Technical secondary school", "College", "Bachelor", expression("">="Master")))  +
  coord_cartesian(ylim = c(13, 16))

ggarrange(p1, p2, p3, p4, 
          p8, p9, p10, p11,
          ncol=4, nrow=2, align = "h")



data.tac <- weekday %>% group_by(ID) %>% 
  dplyr::summarise(tac.mean = mean(tac, na.rm = T), 
                   grade = unique(grade), 
                   gender = unique(gender),
                   income = unique(income), 
                   edu = unique(edu),
                   acips.mean = mean(acips, na.rm = T),
                   erq_cr.mean = mean(erq_cr, na.rm = T), 
                   erq_es.mean = mean(erq_es, na.rm = T),
                   sf12v2_mcs.mean = mean(sf12v2_mcs, na.rm = T),
                   sf12v2_pcs.mean = mean(sf12v2_pcs, na.rm = T),
                   dep.cat = unique(dep_cat, na.rm = T)) %>% 
  group_by(gender, dep.cat, grade) %>%
  dplyr::summarise(tac.mean = mean(tac.mean))

p1 <- ggplot(data.tac, aes(x=grade, y=log(tac.mean), 
                           group=interaction(gender,dep.cat),
                           color=interaction(gender,dep.cat))) + 
  geom_point(size=1.3) + 
  geom_line(size=1.3) +
  theme_bw() + 
  scale_color_discrete(name="Gender and Dep_cat", labels = c("Male NoNDep", "Female NonDep", "Male Dep", "Female Dep")) +
  labs(y = "log(Total Activity Count)", x = "Grade") + 
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14),
        legend.text = element_text(size=12), legend.title = element_text(size=14)) + 
  scale_x_continuous(breaks=1:7, labels=6:12) + 
  coord_cartesian(ylim = c(14, 15))

data.tac <- weekend %>% group_by(ID) %>% 
  dplyr::summarise(tac.mean = mean(tac, na.rm = T), 
                   grade = unique(grade), 
                   gender = unique(gender),
                   income = unique(income), 
                   edu = unique(edu),
                   acips.mean = mean(acips, na.rm = T),
                   erq_cr.mean = mean(erq_cr, na.rm = T), 
                   erq_es.mean = mean(erq_es, na.rm = T),
                   sf12v2_mcs.mean = mean(sf12v2_mcs, na.rm = T),
                   sf12v2_pcs.mean = mean(sf12v2_pcs, na.rm = T),
                   dep.cat = unique(dep_cat, na.rm = T)) %>% 
  group_by(gender, dep.cat, grade) %>%
  dplyr::summarise(tac.mean = mean(tac.mean))

p2 <- ggplot(data.tac, aes(x=grade, y=log(tac.mean), 
                           group=interaction(gender,dep.cat),
                           color=interaction(gender,dep.cat))) + 
  geom_point(size=1.3) + 
  geom_line(size=1.3) +
  theme_bw() + 
  scale_color_discrete(name="Gender and Dep_cat", labels = c("Male NoNDep", "Female NonDep", "Male Dep", "Female Dep")) +
  labs(y = "log(Total Activity Count)", x = "Grade") + 
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14), 
        legend.text = element_text(size=12), legend.title = element_text(size=14)) + 
  scale_x_continuous(breaks=1:7, labels=6:12) +
  coord_cartesian(ylim = c(14, 15))

ggarrange(p1, p2,
          ncol=2, nrow=1, align = "v", common.legend = T, legend = "right")




#######Jiangxi##########

acc.all <- all
weekday <- acc.all[acc.all$days_week != "Saturday" & acc.all$days_week != "Sunday",]
weekend <- acc.all[acc.all$days_week == "Saturday" | acc.all$days_week == "Sunday",]

if (sum(weekday$urban == "urban") > 0){
  data.tac <- weekday %>% group_by(ID) %>% 
    dplyr::summarise(tac.mean = mean(tac, na.rm = T), 
                     grade = unique(grade), 
                     gender = unique(gender),
                     income = unique(income), 
                     edu = unique(edu),
                     acips.mean = mean(acips, na.rm = T),
                     erq_cr.mean = mean(erq_cr, na.rm = T), 
                     erq_es.mean = mean(erq_es, na.rm = T),
                     sf12v2_mcs.mean = mean(sf12v2_mcs, na.rm = T),
                     sf12v2_pcs.mean = mean(sf12v2_pcs, na.rm = T),
                     dep.mean = mean(depression, na.rm = T), 
                     school = unique(school), 
                     judu = unique(judu), 
                     mother.leave = unique(mother.leave), 
                     urban = unique(urban))
} else {
  data.tac <- weekday %>% group_by(ID) %>% 
    dplyr::summarise(tac.mean = mean(tac, na.rm = T), 
                     grade = unique(grade), 
                     gender = unique(gender),
                     income = unique(income), 
                     edu = unique(edu),
                     acips.mean = mean(acips, na.rm = T),
                     erq_cr.mean = mean(erq_cr, na.rm = T), 
                     erq_es.mean = mean(erq_es, na.rm = T),
                     sf12v2_mcs.mean = mean(sf12v2_mcs, na.rm = T),
                     sf12v2_pcs.mean = mean(sf12v2_pcs, na.rm = T),
                     dep.mean = mean(depression, na.rm = T), 
                     school = unique(school))
}

give.n <- function(x){
  return(data.frame(y = quantile(x + 0.3, probs = 0.99), label = length(x)))
  #label = paste0("n = ",length(x))))
}


p1 <- ggplot(data.tac, aes(x=as.factor(gender), y=log(tac.mean), fill = as.factor(gender))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA, width = 0.5) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 4) +
  scale_fill_discrete(name = "Sex", labels = c("Male", "Female")) + 
  #guides(fill=guide_legend(title="Sex")) + 
  labs(y = "log(Total Activity Count) \n Weekday", title = "Gender") + 
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12),  axis.text.x = element_text(size=12), 
        legend.position="none") + 
  scale_x_discrete(breaks=1:2, labels=c("Male", "Female")) + 
  coord_cartesian(ylim = c(13, 16))

p2 <- ggplot(data.tac, aes(x=factor(grade), y=log(tac.mean), fill = as.factor(grade))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 3.5) +
  scale_fill_discrete(name = "Grade", labels = 7:12) + 
  #guides(fill=guide_legend(title="Grade")) + 
  labs(title = "Grade") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), 
        legend.position="none") + 
  scale_x_discrete(breaks=1:7, labels=6:12) + 
  coord_cartesian(ylim = c(13, 16))

p3 <- ggplot(subset(data.tac, !is.na(income)), aes(x=factor(income),  y=log(tac.mean), fill = as.factor(income))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 3.5) +
  scale_fill_discrete(name = "Income", labels = 1:8) + 
  #guides(fill=guide_legend(title="Income")) + 
  labs(title = "Income (¥)") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=10, angle = 30, hjust=1), 
        legend.position="none") +  
  scale_x_discrete(breaks=1:7, labels=c(expression(""<="10k"), "10-30k", "30-50k", "50-100k", 
                                        "100-150k", "150-300k", expression("">="300k"))) + 
  coord_cartesian(ylim = c(13, 16))

p4 <- ggplot(subset(data.tac, !is.na(edu)), aes(x=factor(edu), y=log(tac.mean), fill = as.factor(edu))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 3.5) +
  scale_fill_discrete(name = "Education", labels = 1:9) + 
  #guides(fill=guide_legend(title="Education")) +
  labs(title = "Education") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=10, angle = 30, hjust=1), 
        legend.position="none") +  
  scale_x_discrete(breaks=2:8, labels=c(expression(""<"Primary school"), "Primary school", "Middle school", "High school", 
                                        "Technical secondary school", "College", expression("">="Bachelor")))  +
  coord_cartesian(ylim = c(13, 16))


p5 <- ggplot(subset(data.tac, !is.na(mother.leave)), 
             aes(x=factor(mother.leave), y=log(tac.mean), fill = as.factor(mother.leave))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA, width = 0.5) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 4) + 
  labs(title = "Mother Left Behind") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), 
        legend.position="none") + 
  coord_cartesian(ylim = c(13, 16))

p6 <- ggplot(subset(data.tac, !is.na(judu)), 
             aes(x=factor(judu), y=log(tac.mean), fill = as.factor(judu))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA, width = 0.5) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 4) +
  labs(title = "Lodge") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), 
        legend.position="none") + 
  scale_x_discrete(labels=c("Extern", "Lodge")) + 
  coord_cartesian(ylim = c(13, 16))

p7 <- ggplot(subset(data.tac, !is.na(urban)), 
             aes(x=factor(urban), y=log(tac.mean), fill = as.factor(urban))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA, width = 0.5) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 4) +
  labs(title = "Urban") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), 
        legend.position="none") + 
  scale_x_discrete(labels=c("Rural", "Urban")) + 
  coord_cartesian(ylim = c(13, 16))


ggarrange(p1, p2, p3, p4, p5, p6, p7,
          ncol=4, nrow=2, align = "hv")


if (sum(weekend$urban == "urban") > 0){
  data.tac <- weekend %>% group_by(ID) %>% 
    dplyr::summarise(tac.mean = mean(tac, na.rm = T), 
                     grade = unique(grade), 
                     gender = unique(gender),
                     income = unique(income), 
                     edu = unique(edu),
                     acips.mean = mean(acips, na.rm = T),
                     erq_cr.mean = mean(erq_cr, na.rm = T), 
                     erq_es.mean = mean(erq_es, na.rm = T),
                     sf12v2_mcs.mean = mean(sf12v2_mcs, na.rm = T),
                     sf12v2_pcs.mean = mean(sf12v2_pcs, na.rm = T),
                     dep.mean = mean(depression, na.rm = T), 
                     school = unique(school), 
                     judu = unique(judu), 
                     mother.leave = unique(mother.leave), 
                     urban = unique(urban))
} else {
  data.tac <- weekend %>% group_by(ID) %>% 
    dplyr::summarise(tac.mean = mean(tac, na.rm = T), 
                     grade = unique(grade), 
                     gender = unique(gender),
                     income = unique(income), 
                     edu = unique(edu),
                     acips.mean = mean(acips, na.rm = T),
                     erq_cr.mean = mean(erq_cr, na.rm = T), 
                     erq_es.mean = mean(erq_es, na.rm = T),
                     sf12v2_mcs.mean = mean(sf12v2_mcs, na.rm = T),
                     sf12v2_pcs.mean = mean(sf12v2_pcs, na.rm = T),
                     dep.mean = mean(depression, na.rm = T), 
                     school = unique(school))
}


p8 <- ggplot(data.tac, aes(x=as.factor(gender), y=log(tac.mean), fill = as.factor(gender))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA, width = 0.5) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 4) +
  scale_fill_discrete(name = "Sex", labels = c("Male", "Female")) + 
  #guides(fill=guide_legend(title="Sex")) + 
  labs(y = "log(Total Activity Count) \n Weekend", title = "Gender") + 
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12),  axis.text.x = element_text(size=12), 
        legend.position="none") + 
  scale_x_discrete(breaks=1:2, labels=c("Male", "Female")) + 
  coord_cartesian(ylim = c(13, 16))

p9 <- ggplot(data.tac, aes(x=factor(grade), y=log(tac.mean), fill = as.factor(grade))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 3.5) +
  scale_fill_discrete(name = "Grade", labels = 7:12) + 
  #guides(fill=guide_legend(title="Grade")) + 
  labs(title = "Grade") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), 
        legend.position="none") + 
  scale_x_discrete(breaks=1:7, labels=6:12) + 
  coord_cartesian(ylim = c(13, 16))

p10 <- ggplot(subset(data.tac, !is.na(income)), aes(x=factor(income),  y=log(tac.mean), fill = as.factor(income))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 3.5) +
  scale_fill_discrete(name = "Income", labels = 1:8) + 
  #guides(fill=guide_legend(title="Income")) + 
  labs(title = "Income (¥)") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=10, angle = 30, hjust=1), 
        legend.position="none") +  
  scale_x_discrete(breaks=1:7, labels=c(expression(""<="10k"), "10-30k", "30-50k", "50-100k", 
                                        "100-150k", "150-300k", expression("">="300k"))) + 
  coord_cartesian(ylim = c(13, 16))

p11 <- ggplot(subset(data.tac, !is.na(edu)), aes(x=factor(edu), y=log(tac.mean), fill = as.factor(edu))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 3.5) +
  scale_fill_discrete(name = "Education", labels = 1:9) + 
  #guides(fill=guide_legend(title="Education")) +
  labs(title = "Education") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=10, angle = 30, hjust=1), 
        legend.position="none") +  
  scale_x_discrete(breaks=2:8, labels=c(expression(""<"Primary school"), "Primary school", "Middle school", "High school", 
                                        "Technical secondary school", "College", expression("">="Bachelor")))  +
  coord_cartesian(ylim = c(13, 16))


p12 <- ggplot(subset(data.tac, !is.na(mother.leave)), 
             aes(x=factor(mother.leave), y=log(tac.mean), fill = as.factor(mother.leave))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA, width = 0.5) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 4) + 
  labs(title = "Mother Left Behind") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), 
        legend.position="none") + 
  coord_cartesian(ylim = c(13, 16))

p13 <- ggplot(subset(data.tac, !is.na(judu)), 
             aes(x=factor(judu), y=log(tac.mean), fill = as.factor(judu))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA, width = 0.5) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 4) +
  labs(title = "Lodge") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), 
        legend.position="none") + 
  scale_x_discrete(labels=c("Extern", "Lodge")) + 
  coord_cartesian(ylim = c(13, 16))

p14 <- ggplot(subset(data.tac, !is.na(urban)), 
             aes(x=factor(urban), y=log(tac.mean), fill = as.factor(urban))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA, width = 0.5) + 
  theme_bw() +
  stat_summary(fun.data = give.n, geom = "text",  size = 4) +
  labs(title = "Urban") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), 
        legend.position="none") + 
  scale_x_discrete(labels=c("Rural", "Urban")) + 
  coord_cartesian(ylim = c(13, 16))


ggarrange(p8, p9, p10, p11, 
          p12, p13, p14,
          ncol=4, nrow=2, align = "hv")



data.tac <- weekday %>% group_by(ID) %>% 
  dplyr::summarise(tac.mean = mean(tac, na.rm = T), 
                   grade = unique(grade), 
                   gender = unique(gender),
                   income = unique(income), 
                   edu = unique(edu),
                   acips.mean = mean(acips, na.rm = T),
                   erq_cr.mean = mean(erq_cr, na.rm = T), 
                   erq_es.mean = mean(erq_es, na.rm = T),
                   sf12v2_mcs.mean = mean(sf12v2_mcs, na.rm = T),
                   sf12v2_pcs.mean = mean(sf12v2_pcs, na.rm = T),
                   dep.cat = unique(dep_cat, na.rm = T)) %>% 
  group_by(gender, dep.cat, grade) %>%
  dplyr::summarise(tac.mean = mean(tac.mean))

p1 <- ggplot(data.tac, aes(x=grade, y=log(tac.mean), 
                           group=interaction(gender,dep.cat),
                           color=interaction(gender,dep.cat))) + 
  geom_point(size=1.3) + 
  geom_line(size=1.3) +
  theme_bw() + 
  scale_color_discrete(name="Gender and Dep_cat", labels = c("Male NoNDep", "Female NonDep", "Male Dep", "Female Dep")) +
  labs(y = "log(Total Activity Count)", x = "Grade") + 
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14),
        legend.text = element_text(size=12), legend.title = element_text(size=14)) + 
  scale_x_continuous(breaks=1:7, labels=6:12) + 
  coord_cartesian(ylim = c(14, 15))

data.tac <- weekend %>% group_by(ID) %>% 
  dplyr::summarise(tac.mean = mean(tac, na.rm = T), 
                   grade = unique(grade), 
                   gender = unique(gender),
                   income = unique(income), 
                   edu = unique(edu),
                   acips.mean = mean(acips, na.rm = T),
                   erq_cr.mean = mean(erq_cr, na.rm = T), 
                   erq_es.mean = mean(erq_es, na.rm = T),
                   sf12v2_mcs.mean = mean(sf12v2_mcs, na.rm = T),
                   sf12v2_pcs.mean = mean(sf12v2_pcs, na.rm = T),
                   dep.cat = unique(dep_cat, na.rm = T)) %>% 
  group_by(gender, dep.cat, grade) %>%
  dplyr::summarise(tac.mean = mean(tac.mean))

p2 <- ggplot(data.tac, aes(x=grade, y=log(tac.mean), 
                           group=interaction(gender,dep.cat),
                           color=interaction(gender,dep.cat))) + 
  geom_point(size=1.3) + 
  geom_line(size=1.3) +
  theme_bw() + 
  scale_color_discrete(name="Gender and Dep_cat", labels = c("Male NoNDep", "Female NonDep", "Male Dep", "Female Dep")) +
  labs(y = "log(Total Activity Count)", x = "Grade") + 
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14), 
        legend.text = element_text(size=12), legend.title = element_text(size=14)) + 
  scale_x_continuous(breaks=1:7, labels=6:12) +
  coord_cartesian(ylim = c(14, 15))

ggarrange(p1, p2,
          ncol=2, nrow=1, align = "v", common.legend = T, legend = "right")



data.tac <- weekday %>% group_by(ID) %>% 
  dplyr::summarise(tac.mean = mean(tac, na.rm = T),
                   depression.mean = mean(depression, na.rm = T),
                   grade = unique(grade), 
                   gender = unique(gender),
                   income = unique(income), 
                   edu = unique(edu),
                   mother.leave = unique(mother.leave)) 
#group_by(edu, mother.leave) %>%
#dplyr::summarise(tac.mean = mean(tac.mean), depression.mean = mean(depression.mean))



p1 <- ggplot(data.tac, aes(x=factor(edu), y=depression.mean, 
                           fill=factor(mother.leave))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw()+
  scale_fill_discrete(name = "Mother Left Behind", labels = 0:1) + 
  guides(fill=guide_legend(title="Mother Left Behind")) +  labs(x = "Education", y = "Depression") +
  theme(axis.title.y = element_text(size=16), axis.title.x = element_text(size=14), 
        axis.text.y = element_text(size=14), axis.text.x = element_text(size=12, angle = 30, hjust=1), 
        legend.text = element_text(size=12), legend.title = element_text(size=14)) + 
  scale_x_discrete(breaks=2:8, labels=c(expression(""<"Primary school"), "Primary school", "Middle school", "High school", 
                                        "Technical secondary school", "College", expression("">="Bachelor")))  +
  coord_cartesian(ylim = c(-2, 17)) + 
  geom_text(data = data.tac %>% group_by(edu, mother.leave) %>% 
              summarize(Count = n(), 
                        depression.mean= min(depression.mean)-1),
            aes(label = Count), 
            position = position_dodge(0.95), size=4.5, show.legend = FALSE) 


p2 <- ggplot(data.tac, aes(x=factor(income), y=depression.mean, 
                     fill=factor(mother.leave))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() + 
  scale_fill_discrete(name = "Mother Left Behind", labels = 0:1) + 
  guides(fill=guide_legend(title="Mother Left Behind")) +  labs(x = "Income", y = "Depression") +
  theme(axis.title.y = element_text(size=16), axis.title.x = element_text(size=14), 
        axis.text.y = element_text(size=14), axis.text.x = element_text(size=12, angle = 30, hjust=1), 
        legend.text = element_text(size=12), legend.title = element_text(size=14)) + 
  scale_x_discrete(breaks=1:7, labels=c(expression(""<="10k"), "10-30k", "30-50k", "50-100k", 
                                        "100-150k", "150-300k", expression("">="300k"))) +
  coord_cartesian(ylim = c(-2, 17)) + 
  geom_text(data = data.tac %>% group_by(income, mother.leave) %>% 
              summarize(Count = n(), 
                        depression.mean= min(depression.mean)-1),
            aes(label = Count), 
            position = position_dodge(0.95), size=4.5, show.legend = FALSE) 


ggarrange(p1, p2,
          ncol=2, nrow=1, align = "h", common.legend = T, legend = "right")



####################

cor.data <- data.frame(acips=all$acips, erq_cr=all$erq_cr, erq_es=all$erq_es,
                       sf12v2_mcs=all$sf12v2_mcs, sf12v2_pcs=all$sf12v2_pcs, 
                       anxiety=all$anxiety, depression=all$depression, stress=all$stress)

library(corrplot)
corrplot.mixed(cor(cor.data),
               lower = "number", 
               upper = "circle",
               tl.col = "black",
               lower.col = "black")


####################

p1 <- ggplot(sh_out, aes(x=factor(grade_num), y=acips, fill=factor(new_gender))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  scale_fill_discrete(name = "sex", labels = c("male", "female")) + 
  guides(fill=guide_legend(title="Sex")) + labs(y = "Shanghai", title = "ACIPS") + 
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_blank())

p2 <- ggplot(sh_out, aes(x=factor(grade_num), y=erq_cr, fill=factor(new_gender))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  scale_fill_discrete(name = "sex", labels = c("male", "female")) + 
  guides(fill=guide_legend(title="Sex")) + labs(title = "ERQ_CR") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

p3 <- ggplot(sh_out, aes(x=factor(grade_num), y=erq_es, fill=factor(new_gender))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  scale_fill_discrete(name = "sex", labels = c("male", "female")) + 
  guides(fill=guide_legend(title="Sex")) + labs(title = "ERQ_ES") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

p4 <- ggplot(sh_out, aes(x=factor(grade_num), y=SF12v2_mcs, fill=factor(new_gender))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  scale_fill_discrete(name = "sex", labels = c("male", "female")) + 
  guides(fill=guide_legend(title="Sex")) + labs(title = "SF12v2_mcs") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

p5 <- ggplot(sh_out, aes(x=factor(grade_num), y=SF12v2_pcs, fill=factor(new_gender))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  scale_fill_discrete(name = "sex", labels = c("male", "female")) + 
  guides(fill=guide_legend(title="Sex")) + labs(title = "SF12v2_pcs") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

p6 <- ggplot(sh_out, aes(x=factor(grade_num), y=dass_depression, fill=factor(new_gender))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  scale_fill_discrete(name = "sex", labels = c("male", "female")) + 
  guides(fill=guide_legend(title="Sex")) + labs(title ="Depression") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

p7 <- ggplot(jx_out, aes(x=factor(grade_num), y=acips, fill=factor(new_gender))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  scale_fill_discrete(name = "sex", labels = c("male", "female")) + 
  guides(fill=guide_legend(title="Sex")) +labs(x = "Grade", y = "Jiangxi")+ 
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        plot.title = element_blank())

p8 <- ggplot(jx_out, aes(x=factor(grade_num), y=erq_cr, fill=factor(new_gender))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  scale_fill_discrete(name = "sex", labels = c("male", "female")) + 
  guides(fill=guide_legend(title="Sex")) + labs(x = "Grade") +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 14), 
        plot.title = element_blank())

p9 <- ggplot(jx_out, aes(x=factor(grade_num), y=erq_es, fill=factor(new_gender))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  scale_fill_discrete(name = "sex", labels = c("male", "female")) + 
  guides(fill=guide_legend(title="Sex")) + labs(x = "Grade") +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 14), 
        plot.title = element_blank())

p10 <- ggplot(jx_out, aes(x=factor(grade_num), y=SF12v2_mcs, fill=factor(new_gender))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  scale_fill_discrete(name = "sex", labels = c("male", "female")) + 
  guides(fill=guide_legend(title="Sex")) + labs(x = "Grade") +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 14), 
        plot.title = element_blank())

p11 <- ggplot(jx_out, aes(x=factor(grade_num), y=SF12v2_pcs, fill=factor(new_gender))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  scale_fill_discrete(name = "sex", labels = c("male", "female")) + 
  guides(fill=guide_legend(title="Sex")) + labs(x = "Grade") +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 14), 
        plot.title = element_blank())

p12 <- ggplot(jx_out, aes(x=factor(grade_num), y=dass_depression, fill=factor(new_gender))) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  scale_fill_discrete(name = "sex", labels = c("male", "female")) + 
  guides(fill=guide_legend(title="Sex")) + labs(x = "Grade") +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 14), 
        plot.title = element_blank())


ggarrange(p1, p2, p3, p4, p5, p6, 
          p7, p8, p9, p10, p11, p12, 
          ncol=6, nrow=2, common.legend = T, legend="right", align="h")


######################
sh.mental <- data.frame(city=rep("SH", nrow(sh_out)),  sex=sh_out$new_gender, grade = sh_out$grade_num, 
                        acips = sh_out$acips, erq_cr = sh_out$erq_cr, erq_es = sh_out$erq_es, 
                        sf12v2_mcs = sh_out$SF12v2_mcs, sf12v2_pcs = sh_out$SF12v2_pcs, 
                        depression = sh_out$dass_depression, dep_cat = as.numeric(sh_out$dep_cat))
jx.mental <- data.frame(city=rep("JX", nrow(jx_out)),  sex=jx_out$new_gender, grade = jx_out$grade_num, 
                        acips = jx_out$acips, erq_cr = jx_out$erq_cr, erq_es = jx_out$erq_es, 
                        sf12v2_mcs = jx_out$SF12v2_mcs, sf12v2_pcs = jx_out$SF12v2_pcs, 
                        depression = jx_out$dass_depression, dep_cat = as.numeric(jx_out$dep_cat))

mental <- rbind(sh.mental, jx.mental)
mental$sex <- ifelse(mental$sex == 1, "Male", "Female")
mental$dep_cat <- factor(ifelse(mental$dep_cat != 0, 1, 0))


p1 <- ggplot(mental, aes(x=grade, y=acips, fill=city)) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  scale_fill_discrete(name = "City") + 
  guides(fill=guide_legend(title="City")) + labs(y = "ACIPS") +
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_blank(), 
        plot.title = element_blank()) + 
  facet_wrap(~factor(sex))

p2 <- ggplot(mental, aes(x=grade, y=erq_cr, fill=city)) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  scale_fill_discrete(name = "City") + 
  guides(fill=guide_legend(title="City")) + labs(y = "ERQ_CR") +
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_blank(), 
        plot.title = element_blank()) + 
  facet_wrap(~factor(sex))

p3 <- ggplot(mental, aes(x=grade, y=erq_es, fill=city)) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  scale_fill_discrete(name = "City") + 
  guides(fill=guide_legend(title="City")) + labs(y = "ERQ_ES") +
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_blank(), 
        plot.title = element_blank()) + 
  facet_wrap(~factor(sex))

p4 <- ggplot(mental, aes(x=grade, y=sf12v2_mcs, fill=city)) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  scale_fill_discrete(name = "City") + 
  guides(fill=guide_legend(title="City")) + labs(x = "Grade", y = "SF12V2_mcs") +
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), 
        plot.title = element_blank()) + 
  facet_wrap(~factor(sex))

p5 <- ggplot(mental, aes(x=grade, y=sf12v2_pcs, fill=city)) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  scale_fill_discrete(name = "City") + 
  guides(fill=guide_legend(title="City")) + labs(x = "Grade", y = "SF12V2_pcs") +
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), 
        plot.title = element_blank()) + 
  facet_wrap(~factor(sex))


p6 <- ggplot(mental, aes(x=grade, y=depression, fill=city)) + 
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  theme_bw() +
  scale_fill_discrete(name = "City") + 
  guides(fill=guide_legend(title="City")) + labs(x = "Grade", y = "Depression") +
  theme(axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), 
        plot.title = element_blank()) + 
  facet_wrap(~factor(sex))

ggarrange(p1, p2, p3, p4, p5, p6,
          ncol=3, nrow=2, common.legend = T, legend="right", align="h")




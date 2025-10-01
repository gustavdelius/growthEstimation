age_at_length<-readRDS("/Users/jessicawestworth/Desktop/Mizer Work/Projects/mizerEcopathedits/inst/extdata/Celtic_Sea_Size_at_Age_Data.rds")
length_freq<-readRDS("/Users/jessicawestworth/Desktop/Mizer Work/Projects/mizerEcopathedits/inst/extdata/cs_survey.rds")

library(dplyr)
surveys <-length_freq_data %>%
    mutate(Length=length,
           Month=Month/12-(1/24), ##assume middle of the month
           survey_date= Month)%>%
    group_by(species,Scientific_name, Length, survey_date, dl) %>%
    summarise(count=sum(catch, na.rm = TRUE), .groups="drop")
surveys<-na.omit(surveys)
length_freq<-subset(surveys, species=="Cod")

age_at_length<-age_at_length %>%
    mutate(Length=LngtClass,
           Month=Month/12-(1/24),
           survey_date= Month,
           K = Age)%>%
    group_by(Scientific_name, Length, survey_date, K) %>%
    summarise(count=sum(CANoAtLngt, na.rm = TRUE), .groups="drop")
age_at_length<-na.omit(age_at_length)
age_at_length<-subset(age_at_length, Scientific_name=="Gadus morhua")


#Plotting length_freq.R

pars <- list(
    k = 0.2,
    L_inf = 90,
    d = 0.1,
    m = 30,
    r = 0,
    annuli_date=0,
    annuli_min_age = 0,
    spawning_mu = 0.5,
    spawning_kappa = 70,
    l50 = 10,      # initial value for selectivity midpoint
    ratio = 0.5    # initial value for sigmoidal selectivity ratio
)

fit<-fit_tmb_nll (
    pars = pars,
    age_at_length = age_at_length,
    length_freq = length_freq
)

fit$pars
fit$opt$objective
spectra<-getPeriodicNumberDensity(
        pars=fit$pars,
        l_max=(max(length_freq$Length)*1.1),
        Delta_l = 1,
        t_max = 1,
        Delta_t = 1/24)

selectivity<-vector(length=((max(length_freq$Length)*1.1)+1))

for(l in 1:((max(length_freq$Length)*1.1)+1)){
c1 <- 1
sr <- fit$pars$l50 * (c1 - fit$pars$ratio)
s1 <- fit$pars$l50 * log(3.0) / sr
s2 <- s1 / fit$pars$l50;
selectivity[l] <- 1/(c1 + exp(s1 - s2 * l))}

l<-(0:(max(length_freq$Length)*1.1)+1)
l
#plot(selectivity~l, type="l")

model<-sweep(spectra, 2, selectivity, '*')
# Join survey data to this sequence

survey1<-length_freq %>%
    filter(near(survey_date, 1/24))
survey1$normcount<-survey1$count/(sum(survey1$count))
survey1 <- all_lengths %>%
    left_join(survey1 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))

survey2<-length_freq %>%
    filter(near(survey_date, 3/24))
survey2$normcount<-survey2$count/(sum(survey2$count))
survey2 <- all_lengths %>%
    left_join(survey2 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))

survey3<-length_freq %>%
    filter(near(survey_date, 5/24))
survey3$normcount<-survey3$count/(sum(survey3$count))
survey3 <- all_lengths %>%
    left_join(survey3 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))
survey4<-length_freq %>%
    filter(near(survey_date, 7/24))
survey4$normcount<-survey4$count/(sum(survey4$count))
survey4 <- all_lengths %>%
    left_join(survey4 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))
survey5<-length_freq %>%
    filter(near(survey_date, 9/24))
survey5$normcount<-survey5$count/(sum(survey5$count))
survey5 <- all_lengths %>%
    left_join(survey5 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))
survey6<-length_freq %>%
    filter(near(survey_date, 11/24))
survey6$normcount<-survey6$count/(sum(survey6$count))
survey6 <- all_lengths %>%
    left_join(survey6 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))
survey7<-length_freq %>%
    filter(near(survey_date, 13/24))
survey7$normcount<-survey7$count/(sum(survey7$count))
survey7 <- all_lengths %>%
    left_join(survey7 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))
survey8<-length_freq %>%
    filter(near(survey_date, 15/24))
survey8$normcount<-survey8$count/(sum(survey8$count))
survey8 <- all_lengths %>%
    left_join(survey8 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))
survey9<-length_freq %>%
    filter(near(survey_date, 17/24))
survey9$normcount<-survey9$count/(sum(survey9$count))
survey9 <- all_lengths %>%
    left_join(survey9 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))
survey10<-length_freq %>%
    filter(near(survey_date, 19/24))
survey10$normcount<-survey10$count/(sum(survey10$count))
survey10 <- all_lengths %>%
    left_join(survey10 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))
survey11<-length_freq %>%
    filter(near(survey_date, 21/24))
survey11$normcount<-survey11$count/(sum(survey11$count))
survey11 <- all_lengths %>%
    left_join(survey11 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))
survey12<-length_freq %>%
    filter(near(survey_date, 23/24))
survey12$normcount<-survey12$count/(sum(survey12$count))
survey12 <- all_lengths %>%
    left_join(survey12 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))

par(mfrow = c(3, 4))

plot(survey1$normcount~survey1$Length, type="l", main="Survey 1")
points(model[1,]/sum(model[1,])~l, type="l", col="blue")
legend("topleft", legend=c("survey","model"),fill=c("black","blue"), cex =0.5)

plot(survey2$normcount~survey2$Length, type="l", main="Survey 2")
points(model[3,]/sum(model[3,])~l, type="l", col="blue")
legend("topleft", legend=c("survey","model"),fill=c("black","blue"), cex =0.5)

plot(survey3$normcount~survey3$Length, type="l", main="Survey 3")
points(model[5,]/sum(model[5,])~l, type="l", col="blue")
legend("topleft", legend=c("survey","model"),fill=c("black","blue"), cex = 0.5)

plot(survey4$normcount~survey4$Length, type="l", main="Survey 4")
points(model[7,]/sum(model[7,])~l, type="l", col="blue")
legend("topleft", legend=c("survey","model"),fill=c("black","blue"), cex =0.5)

plot(survey5$normcount~survey5$Length, type="l", main="Survey 5")
points(model[9,]/sum(model[9,])~l, type="l", col="blue")
legend("topleft", legend=c("survey","model"),fill=c("black","blue"), cex =0.5)

plot(survey6$normcount~survey6$Length, type="l", main="Survey 6")
points(model[11,]/sum(model[11,])~l, type="l", col="blue")
legend("topleft", legend=c("survey","model"),fill=c("black","blue"), cex =0.5)

plot(survey7$normcount~survey7$Length, type="l", main="Survey 7")
points(model[13,]/sum(model[13,])~l, type="l", col="blue")
legend("topleft", legend=c("survey","model"),fill=c("black","blue"), cex =0.5)

plot(survey8$normcount~survey8$Length, type="l", main="Survey 8")
points(model[15,]/sum(model[15,])~l, type="l", col="blue")
legend("topleft", legend=c("survey","model"),fill=c("black","blue"), cex =0.5)

plot(survey9$normcount~survey9$Length, type="l", main="Survey 9")
points(model[17,]/sum(model[17,])~l, type="l", col="blue")
legend("topleft", legend=c("survey","model"),fill=c("black","blue"), cex =0.5)

plot(survey10$normcount~survey10$Length, type="l", main="Survey 10")
points(model[19,]/sum(model[19,])~l, type="l", col="blue")
legend("topleft", legend=c("survey","model"),fill=c("black","blue"), cex =0.5)

plot(survey11$normcount~survey11$Length, type="l", main="Survey 11")
points(model[21,]/sum(model[21,])~l, type="l", col="blue")
legend("topleft", legend=c("survey","model"),fill=c("black","blue"), cex =0.5)

plot(survey12$normcount~survey12$Length, type="l", main="Survey 12")
points(model[23,]/sum(model[23,])~l, type="l", col="blue")
legend("topleft", legend=c("survey","model"),fill=c("black","blue"), cex =0.5)



#new pars
pars <- list(
    k = 0.5,
    L_inf = 200,
    d = 0.1,
    m = 50,
    spawning_mu = 0.2,
    spawning_kappa = 0.2,
    l50 = 10,      # initial value for selectivity midpoint
    ratio = 0.5    # initial value for sigmoidal selectivity ratio
)

fit<-fit_tmb_length_freq (
    pars = pars,
    length_freq = length_freq)
fit$pars
fit$opt$objective
spectra<-getPeriodicNumberDensity(
    pars=fit$pars,
    l_max=(max(length_freq$Length)*1.1),
    Delta_l = 1,
    t_max = 10,
    Delta_t = 0.05)

selectivity<-vector(length=((max(length_freq$Length)*1.1)+1))

for(l in 1:((max(length_freq$Length)*1.1)+1)){
    c1 <- 1
    sr <- fit$pars$l50 * (c1 - fit$pars$ratio)
    s1 <- fit$pars$l50 * log(3.0) / sr
    s2 <- s1 / fit$pars$l50;
    selectivity[l] <- 1/(c1 + exp(s1 - s2 * l))}

l<-(0:(max(length_freq$Length)*1.1)+1)
l
#plot(selectivity~l, type="l")

model<-sweep(spectra, 2, selectivity, '*')
# Join survey data to this sequence


survey_all<-length_freq
survey_all$normcount<-survey_all$count/(sum(survey_all$count))
all_lengths <- tibble(Length = (0:(max(length_freq$Length)*1.1)+1))
survey_all <- all_lengths %>%
    left_join(survey_all %>%
                  group_by(Length) %>%
                  summarise(normcount = sum(normcount)), by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))


survey_all<-length_freq
survey_all$normcount<-survey_all$count/(sum(survey_all$count))
all_lengths <- tibble(Length = (0:(max(length_freq$Length)*1.1)+1))
survey_all <- all_lengths %>%
    left_join(survey_all %>%
                  group_by(Length) %>%
                  summarise(normcount = sum(normcount)), by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))

survey1<-length_freq %>%
    filter(near(survey_date, 1/24))
survey1$normcount<-survey1$count/(sum(survey1$count))
survey1 <- all_lengths %>%
    left_join(survey1 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))

survey2<-length_freq %>%
    filter(near(survey_date, 3/24))
survey2$normcount<-survey2$count/(sum(survey2$count))
survey2 <- all_lengths %>%
    left_join(survey2 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))

survey3<-length_freq %>%
    filter(near(survey_date, 5/24))
survey3$normcount<-survey3$count/(sum(survey3$count))
survey3 <- all_lengths %>%
    left_join(survey3 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))
survey4<-length_freq %>%
    filter(near(survey_date, 7/24))
survey4$normcount<-survey4$count/(sum(survey4$count))
survey4 <- all_lengths %>%
    left_join(survey4 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))
survey5<-length_freq %>%
    filter(near(survey_date, 9/24))
survey5$normcount<-survey5$count/(sum(survey5$count))
survey5 <- all_lengths %>%
    left_join(survey5 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))
survey6<-length_freq %>%
    filter(near(survey_date, 11/24))
survey6$normcount<-survey6$count/(sum(survey6$count))
survey6 <- all_lengths %>%
    left_join(survey6 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))
survey7<-length_freq %>%
    filter(near(survey_date, 13/24))
survey7$normcount<-survey7$count/(sum(survey7$count))
survey7 <- all_lengths %>%
    left_join(survey7 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))
survey8<-length_freq %>%
    filter(near(survey_date, 15/24))
survey8$normcount<-survey8$count/(sum(survey8$count))
survey8 <- all_lengths %>%
    left_join(survey8 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))
survey9<-length_freq %>%
    filter(near(survey_date, 17/24))
survey9$normcount<-survey9$count/(sum(survey9$count))
survey9 <- all_lengths %>%
    left_join(survey9 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))
survey10<-length_freq %>%
    filter(near(survey_date, 19/24))
survey10$normcount<-survey10$count/(sum(survey10$count))
survey10 <- all_lengths %>%
    left_join(survey10 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))
survey11<-length_freq %>%
    filter(near(survey_date, 21/24))
survey11$normcount<-survey11$count/(sum(survey11$count))
survey11 <- all_lengths %>%
    left_join(survey11 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))
survey12<-length_freq %>%
    filter(near(survey_date, 23/24))
survey12$normcount<-survey12$count/(sum(survey12$count))
survey12 <- all_lengths %>%
    left_join(survey12 %>% group_by(Length) %>% summarise(normcount = sum(normcount)),
              by = "Length") %>%
    mutate(normcount = ifelse(is.na(normcount), 0, normcount))

par(mfrow = c(3, 4))

plot(survey1$normcount~survey1$Length, type="l", main="Survey 1")
points(model[1,]/sum(model[1,])~l, type="l", col="red")
legend("topleft", legend=c("survey","model"),fill=c("black","red"), cex =0.5)

plot(survey2$normcount~survey2$Length, type="l", main="Survey 2")
points(model[3,]/sum(model[3,])~l, type="l", col="red")
legend("topleft", legend=c("survey","model"),fill=c("black","red"), cex =0.5)

plot(survey3$normcount~survey3$Length, type="l", main="Survey 3")
points(model[5,]/sum(model[5,])~l, type="l", col="red")
legend("topleft", legend=c("survey","model"),fill=c("black","red"), cex = 0.5)

plot(survey4$normcount~survey4$Length, type="l", main="Survey 4")
points(model[7,]/sum(model[7,])~l, type="l", col="red")
legend("topleft", legend=c("survey","model"),fill=c("black","red"), cex =0.5)

plot(survey5$normcount~survey5$Length, type="l", main="Survey 5")
points(model[9,]/sum(model[9,])~l, type="l", col="red")
legend("topleft", legend=c("survey","model"),fill=c("black","red"), cex =0.5)

plot(survey6$normcount~survey6$Length, type="l", main="Survey 6")
points(model[11,]/sum(model[11,])~l, type="l", col="red")
legend("topleft", legend=c("survey","model"),fill=c("black","red"), cex =0.5)

plot(survey7$normcount~survey7$Length, type="l", main="Survey 7")
points(model[13,]/sum(model[13,])~l, type="l", col="red")
legend("topleft", legend=c("survey","model"),fill=c("black","red"), cex =0.5)

plot(survey8$normcount~survey8$Length, type="l", main="Survey 8")
points(model[15,]/sum(model[15,])~l, type="l", col="red")
legend("topleft", legend=c("survey","model"),fill=c("black","red"), cex =0.5)

plot(survey9$normcount~survey9$Length, type="l", main="Survey 9")
points(model[17,]/sum(model[17,])~l, type="l", col="red")
legend("topleft", legend=c("survey","model"),fill=c("black","red"), cex =0.5)

plot(survey10$normcount~survey10$Length, type="l", main="Survey 10")
points(model[19,]/sum(model[19,])~l, type="l", col="red")
legend("topleft", legend=c("survey","model"),fill=c("black","red"), cex =0.5)

plot(survey11$normcount~survey11$Length, type="l", main="Survey 11")
points(model[21,]/sum(model[21,])~l, type="l", col="red")
legend("topleft", legend=c("survey","model"),fill=c("black","red"), cex =0.5)

plot(survey12$normcount~survey12$Length, type="l", main="Survey 12")
points(model[23,]/sum(model[23,])~l, type="l", col="red")
legend("topleft", legend=c("survey","model"),fill=c("black","red"), cex =0.5)



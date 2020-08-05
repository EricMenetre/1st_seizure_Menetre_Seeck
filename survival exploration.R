# Survival analysis
library(tidyverse)
library(lubridate)
library(survival)
library(survminer)

#################################################################
##    Data exploration by survival plots and log rank tests    ##
#################################################################


# data importation and transformation
#' delay was calculated, missing values were excluded, the unknown diagnosis were classified
#' as censored, diagnoses were dichotomize in epilepsy, not epilepsy, only the positive delays
#' were taken into account, relapses were dichotomize in relapses or no relapse. For epileptic
#' patients only relapses under treatment were considered as relapses. Data were transformed
#' in factor and in numerical data. Only delays under 2 years (730 days) were used.
data_surv <- data_full_cleaned%>%
  select(code.pat, date.1st.cs.su, date.dx.final, dx.event, n.tot.recidives, motif_relapse, detail_relapse)%>%
  filter(!is.na(date.1st.cs.su), !is.na(date.dx.final))%>%
  mutate(date.1st.cs.su = as_date(date.1st.cs.su),
         date.dx.final  = as_date(date.dx.final),
         censored = ifelse(is.na(dx.event) | dx.event == "Unknown", "censored", "diagnosed"),
         dx.dicho = ifelse(dx.event == "Epilepsy", "epilepsy", "other"),
         delay_dx_final = date.dx.final-date.1st.cs.su)%>%
  filter(delay_dx_final >= 0)%>%
  select(-date.dx.final, -date.1st.cs.su, -dx.event)%>%
  mutate(detail_relapse = ifelse(is.na(detail_relapse), "other", detail_relapse),
         relapses = case_when(n.tot.recidives == 0 ~ "no relapse",
                              detail_relapse == "ttt ok et bonne compliance" | detail_relapse == "pas de ttt" | detail_relapse == "ttt ok taux infra ou ttt inadequat" ~ "relapse",
                              dx.dicho == "epilepsy" & n.tot.recidives > 0 & detail_relapse != "ttt ok et bonne compliance" ~ "no relapse",
                              TRUE ~ "relapse"))%>%
  select(-n.tot.recidives, -motif_relapse, -detail_relapse)%>%
  mutate(censored = as.factor(censored),
         dx.dicho = as.factor(dx.dicho),
         relapses = as.factor(relapses),
         delay_dx_final = as.numeric(delay_dx_final))%>%
  mutate(censored = as.numeric(censored),
         dx.dicho = as.numeric(dx.dicho),
         relapses = as.numeric(relapses))%>%
  filter(delay_dx_final <= 730)

# Fitting of the survival model
#' This model investigates the differences in time to reach the diagnosis among epileptic
#' and non-epileptic patients.
fit <- survfit(Surv(delay_dx_final, censored) ~ dx.dicho, data = data_surv)
print(fit)
summary(fit)$table

#' Graphically, it seems that non-epileptic patients tend to reach their diagnosis later
#' than epileptic patients --> Probably because the patients comes to the ER, where a 
#' first seizure is suspected and it takes time to rule out epilepsy. Another explanation
#' could be that since the epilepsy-like event are less life threatning, it is less urgent
#' to investigate (often syncopes).
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  pval = TRUE,
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in days",   # customize X axis label.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("Epilepsy", "Not epilepsy"),    # change legend labels.
  palette = 
    c("#E7B800", "#2E9FDF"), # custom color palettes.
  fun = "event"
  )
surv_diff <- survdiff(Surv(delay_dx_final, censored) ~ dx.dicho + relapses, data = data_surv)
surv_diff

#' By adding the relapse information, it appears that the patients who relapse tend to recieve
#' their diagnosis later than those who don't. The effect mignt be explained by the difficulty
#' to diagnose some patients. Patients who are hard to diagnose have more relapses until we
#' are sure about the diagnosis and start a treatment. To see the real effect of the delay, 
#' it might be more appropriate to use another variable, such as the delay until the 
#' first consultation.
fit_2 <- survfit(Surv(delay_dx_final, censored) ~ dx.dicho + relapses, data = data_surv)

ggsurv <- ggsurvplot(fit_2, fun = "event", conf.int = TRUE,
                     ggtheme = theme_bw(),
                     legend.labs = c("epi - sz free", "epi - relapse", "non-epi - no relapse", "non-epi - relapse"))
ggsurv$plot +theme_bw() + 
  theme (legend.position = "right")+
  facet_grid(dx.dicho ~ 1)


data_surv_cs <- data_full_cleaned%>%
  select(code.pat, date.1st.cs.su, date.dx.final, dx.event, n.tot.recidives, motif_relapse, detail_relapse,
         date.1st.cs.epi)%>%
  filter(!is.na(date.1st.cs.su), !is.na(date.1st.cs.epi))%>%
  mutate(date.1st.cs.su = as_date(date.1st.cs.su),
         date.dx.final  = as_date(date.dx.final),
         date.1st.cs.epi = as_date(date.1st.cs.epi),
         censored = ifelse(is.na(dx.event) | dx.event == "Unknown", "censored", "diagnosed"),
         dx.dicho = ifelse(dx.event == "Epilepsy", "epilepsy", "other"),
         delay_dx_final = date.dx.final-date.1st.cs.su,
         delay_cs_epi = date.1st.cs.epi - date.1st.cs.su)%>%
  filter(delay_cs_epi >= 0)%>%
  select(-date.dx.final, -date.1st.cs.su, -dx.event, -date.1st.cs.epi)%>%
  mutate(detail_relapse = ifelse(is.na(detail_relapse), "other", detail_relapse),
         relapses = case_when(n.tot.recidives == 0 ~ "no relapse",
                              detail_relapse == "ttt ok et bonne compliance" | detail_relapse == "pas de ttt" | detail_relapse == "ttt ok taux infra ou ttt inadequat" ~ "relapse",
                              dx.dicho == "epilepsy" & n.tot.recidives > 0 & detail_relapse != "ttt ok et bonne compliance" ~ "no relapse",
                              TRUE ~ "relapse"))%>%
  select(-n.tot.recidives, -motif_relapse, -detail_relapse)%>%
  mutate(censored = as.factor(censored),
         dx.dicho = as.factor(dx.dicho),
         relapses = as.factor(relapses),
         delay_dx_final = as.numeric(delay_dx_final))%>%
  mutate(censored = as.numeric(censored),
         dx.dicho = as.numeric(dx.dicho),
         relapses = as.numeric(relapses))%>%
  filter(delay_cs_epi <= 730)


#' It seems that non-epileptic patients are seen in priority compared to epileptic patients.
#' this counter intuitive result might be explained by the fact that easy epileptic cases
#' are seen and treated in the ED, implying that the follow-up is not so urgent.
fit_cs <- survfit(Surv(delay_cs_epi, censored) ~ dx.dicho, data = data_surv_cs)
print(fit_cs)
summary(fit_cs)$table

ggsurvplot(
  fit_cs,                     # survfit object with calculated statistics.
  pval = TRUE,
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in days",   # customize X axis label.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("Epilepsy", "Not epilepsy"),    # change legend labels.
  palette = 
    c("#E7B800", "#2E9FDF"), # custom color palettes.
  fun = "event"
)
surv_diff <- survdiff(Surv(delay_dx_final, censored) ~ dx.dicho + relapses, data = data_surv)
surv_diff

data_full_cleaned%>%
  mutate(delay_cs = as_date(date.1st.cs.epi)-as_date(date.1st.cs.su),
         dx.dicho = ifelse(dx.event != "Epilepsy" | is.na(dx.event), "non-epi", "epi"))%>%
  filter(delay_cs >= 0, delay_cs < 60)%>%
  group_by(dx.dicho)%>%
  summarise(mean_delay = mean(delay_cs),
            median_delay = median(delay_cs))

#' There does not seems to be an effect according to the delay for the epileptic patients
#' but more for the non-epileptic patients.
fit_cs_2 <- survfit(Surv(delay_cs_epi, censored) ~ dx.dicho+relapses, data = data_surv_cs)
print(fit_cs_2)
summary(fit)$table

ggsurv <- ggsurvplot(fit_cs_2, fun = "event", conf.int = TRUE,
                     pval = TRUE,
                     ggtheme = theme_bw(),
                     legend.labs = c("epi - sz free", "epi - relapse", "non-epi - no relapse", "non-epi - relapse"))
ggsurv$plot +theme_bw() + 
  theme (legend.position = "right")+
  facet_grid(dx.dicho ~ 1)

#' No effect, probably because patients are treated in the ED. Let's try with
#' treatment delay.
data_surv_cs_epi <- data_surv_cs%>%filter(dx.dicho == 1)
fit_cs_3 <- survfit(Surv(delay_cs_epi, censored) ~ relapses, data = data_surv_cs_epi)
print(fit_cs_3)

ggsurvplot(fit_cs_3, fun = "event",
                     pval = TRUE,
                     legend.labs = c("no relapse", "relapse"),
           conf.int = TRUE,         # show confidence intervals for 
           # point estimaes of survival curves.
           conf.int.style = "step",  # customize style of confidence intervals
           xlab = "Time in days",   # customize X axis label.
           ggtheme = theme_light(), # customize plot and risk table with a theme.
           # in legend of risk table.
           ncensor.plot = TRUE,      # plot the number of censored subjects at time t
           surv.median.line = "hv",  # add the median survival pointer.
           )

#' The delay between the first seizure and the instauration of treatment does not seems
#' to impact the presence of relapses.
data_surv_ttt <- data_full_cleaned%>%
  select(code.pat, date.1st.cs.su, date.dx.final, dx.event, n.tot.recidives, motif_relapse, detail_relapse,
         date.1st.cs.epi, date.debut.tt)%>%
  filter(!is.na(date.1st.cs.su),!is.na(date.debut.tt),
         dx.event == "Epilepsy")%>%
  mutate(date.1st.cs.su = as_date(date.1st.cs.su),
         date.dx.final  = as_date(date.dx.final),
         date.1st.cs.epi = as_date(date.1st.cs.epi),
         date.debut.tt = as_date(date.debut.tt),
         censored = ifelse(is.na(dx.event) | dx.event == "Unknown", "censored", "diagnosed"),
         dx.dicho = ifelse(dx.event == "Epilepsy", "epilepsy", "other"),
         delay_dx_final = date.dx.final-date.1st.cs.su,
         delay_cs_epi = date.1st.cs.epi - date.1st.cs.su,
         delay_tt = date.debut.tt - date.1st.cs.su)%>%
  filter(delay_tt >= 0)%>%
  select(-date.dx.final, -date.1st.cs.su, -dx.event, -date.1st.cs.epi, -date.debut.tt)%>%
  mutate(detail_relapse = ifelse(is.na(detail_relapse), "other", detail_relapse),
         relapses = case_when(n.tot.recidives == 0 ~ "no relapse",
                              detail_relapse == "ttt ok et bonne compliance" | detail_relapse == "pas de ttt" | detail_relapse == "ttt ok taux infra ou ttt inadequat" ~ "relapse",
                              dx.dicho == "epilepsy" & n.tot.recidives > 0 & detail_relapse != "ttt ok et bonne compliance" ~ "no relapse",
                              TRUE ~ "relapse"))%>%
  select(-n.tot.recidives, -motif_relapse, -detail_relapse)%>%
  mutate(censored = as.factor(censored),
         dx.dicho = as.factor(dx.dicho),
         relapses = as.factor(relapses),
         delay_dx_final = as.numeric(delay_dx_final),
         delay_tt = as.numeric(delay_tt))%>%
  mutate(censored = as.numeric(censored),
         dx.dicho = as.numeric(dx.dicho),
         relapses = as.numeric(relapses))%>%
  filter(delay_tt <=730)


fit_tt <- survfit(Surv(delay_cs_epi, censored) ~ relapses, data = data_surv_ttt)
print(fit_tt)
summary(fit_tt)$table

ggsurvplot(fit_tt, fun = "event", conf.int = TRUE,
                     pval = TRUE,
                     ggtheme = theme_bw(),
                     legend.labs = c("no relapse", "relapse"))

#################################################################
##                 Data analysis by Cox models                 ##
#################################################################

# Set working directory
setwd("~/Documents/Estatistica_mestrado/")

# Load packages
library(openxlsx)
library(dplyr)
library(car)

# Load data
zkv_data <- read.xlsx("metadados_zikv.xlsx")

# Prepare data
vars_factor <- c("Group", "Sexo", "Raça/cor", "Escolaridade", "Renda",
                 "Trim..Sint..de.Zika", "Contato.Álcool", "Contato.Drogas",
                 "Contato.Tabaco/Fumo", "Mãe.vacinada.para.FA.Antes.da.Gest.",
                 "História.Familiar.Consang.", "História.Familiar.com.Malfor.")

for (var in vars_factor) {
  zkv_data[[var]] <- as.factor(zkv_data[[var]])
}

vars_numeric <- c("Idade.Gest..Parto.(Sem.)", "PC.ao.nascer", "Peso.ao.nascer", "Estatura.ao.nascer")

for (var in vars_numeric) {
  zkv_data[[var]] <- as.numeric(zkv_data[[var]])
}


#___ Sociodemographic characteristics between groups ___#
## Sex
assoc_sex <- chisq.test(table(zkv_data$Group, zkv_data$Sexo))
print(assoc_sex)
table(zkv_data$Group, zkv_data$Sexo)

## Skin color
assoc_skin_color <- fisher.test(table(zkv_data$Group, zkv_data$`Raça/cor`), simulate.p.value = TRUE, B = 10000)
print(assoc_skin_color)
table(zkv_data$Group, zkv_data$`Raça/cor`)

## Maternal education level
assoc_education <- fisher.test(table(zkv_data$Group, zkv_data$Escolaridade), simulate.p.value = TRUE, B = 10000)
print(assoc_education)
table(zkv_data$Group, zkv_data$Escolaridade)

## Monthly family income
assoc_income <- fisher.test(table(zkv_data$Group, zkv_data$Renda), simulate.p.value = TRUE, B = 10000)
print(assoc_income)
table(zkv_data$Group, zkv_data$Renda)


#___ Gestational characteristics between groups ___#
## Trimester of infection
assoc_trim <- chisq.test(table(zkv_data$Group, zkv_data$`Trim..Sint..de.Zika`))
print(assoc_trim)
table(zkv_data$Group, zkv_data$`Trim..Sint..de.Zika`)

## Alcohol exposure during pregnancy
assoc_alcohol <- chisq.test(table(zkv_data$Group, zkv_data$`Contato.Álcool`))
print(assoc_alcohol)
table(zkv_data$Group, zkv_data$`Contato.Álcool`)

## Recreational drugs exposure during pregnancy
assoc_drugs <- fisher.test(table(zkv_data$Group, zkv_data$`Contato.Drogas`), simulate.p.value = TRUE, B = 10000)
print(assoc_drugs)
table(zkv_data$Group, zkv_data$`Contato.Drogas`)

## Smoke exposure during pregnancy
assoc_smoke <- fisher.test(table(zkv_data$Group, zkv_data$`Contato.Tabaco/Fumo`), simulate.p.value = TRUE, B = 10000)
print(assoc_smoke)
table(zkv_data$Group, zkv_data$`Contato.Tabaco/Fumo`)

## Gestational age at birth
shapiro.test(zkv_data$`Idade.Gest..Parto.(Sem.)`[zkv_data$Group == "1"])
shapiro.test(zkv_data$`Idade.Gest..Parto.(Sem.)`[zkv_data$Group == "0"])
assoc_gestational_age <- wilcox.test(`Idade.Gest..Parto.(Sem.)` ~ Group, data = zkv_data)
print(assoc_gestational_age)

summary_gestational_age <- aggregate(`Idade.Gest..Parto.(Sem.)` ~ Group, data = zkv_data,
                                     FUN = function(x) c(median=median(x, na.rm=TRUE),
                                                         Q1 = quantile(x, 0.25, na.rm=TRUE),                                 
                                                         Q3 = quantile(x, 0.75, na.rm=TRUE),
                                                         IQR=IQR(x, na.rm=TRUE),
                                                         min=min(x, na.rm=TRUE), max=max(x, na.rm=TRUE)))
print(summary_gestational_age)


#___ Clinical characteristics between groups ___#
## Maternal yellow fever vaccine
filtered_YFv <- subset(zkv_data, `Mãe.vacinada.para.FA.Antes.da.Gest.` %in% c(1, 2))
filtered_YFv$`Mãe.vacinada.para.FA.Antes.da.Gest.` <- droplevels(filtered_YFv$`Mãe.vacinada.para.FA.Antes.da.Gest.`)
assoc_FAv <- chisq.test(table(filtered_YFv$Group, filtered_YFv$`Mãe.vacinada.para.FA.Antes.da.Gest.`))
print(assoc_FAv)
table(filtered_YFv$Group, filtered_YFv$`Mãe.vacinada.para.FA.Antes.da.Gest.`)

## Family history of consanguinity
assoc_hist_cosang <- fisher.test(table(zkv_data$Group, zkv_data$`História.Familiar.Consang.`), simulate.p.value = TRUE, B = 10000)
print(assoc_hist_cosang)
table(zkv_data$Group, zkv_data$`História.Familiar.Consang.`)

## Family history of malformations
filtered_hist_malfor <- subset(zkv_data, `História.Familiar.com.Malfor.` %in% c(1, 2))
filtered_hist_malfor$`História.Familiar.com.Malfor.` <- droplevels(filtered_hist_malfor$`História.Familiar.com.Malfor.`)
assoc_hmalf <- fisher.test(table(filtered_hist_malfor$Group, filtered_hist_malfor$`História.Familiar.com.Malfor.`), simulate.p.value = TRUE, B = 10000)
print(assoc_hmalf)
table(filtered_hist_malfor$Group, filtered_hist_malfor$`História.Familiar.com.Malfor.`)

## Cephalic perimeter
shapiro.test(zkv_data$PC.ao.nascer[zkv_data$Group == "1"])
shapiro.test(zkv_data$PC.ao.nascer[zkv_data$Group == "0"])
assoc_cephalic_perim <- wilcox.test(PC.ao.nascer ~ Group, data = zkv_data)
print(assoc_cephalic_perim)

summary_cephalic_perim <- aggregate(PC.ao.nascer ~ Group, data = zkv_data,
                                    FUN = function(x) c(median=median(x, na.rm=TRUE),
                                                        Q1 = quantile(x, 0.25, na.rm=TRUE),                                 
                                                        Q3 = quantile(x, 0.75, na.rm=TRUE),
                                                        IQR=IQR(x, na.rm=TRUE),
                                                        min=min(x, na.rm=TRUE), max=max(x, na.rm=TRUE)))
print(summary_cephalic_perim)

## Weight
shapiro.test(zkv_data$Peso.ao.nascer[zkv_data$Group == "1"])
shapiro.test(zkv_data$Peso.ao.nascer[zkv_data$Group == "0"])
assoc_weight <- wilcox.test(Peso.ao.nascer ~ Group, data = zkv_data)
print(assoc_weight)

summary_weight <- aggregate(Peso.ao.nascer ~ Group, data = zkv_data,
                            FUN = function(x) c(median=median(x, na.rm=TRUE),
                                                Q1 = quantile(x, 0.25, na.rm=TRUE),                                 
                                                Q3 = quantile(x, 0.75, na.rm=TRUE),
                                                IQR=IQR(x, na.rm=TRUE),
                                                min=min(x, na.rm=TRUE), max=max(x, na.rm=TRUE)))
print(summary_weight)

## Height
shapiro.test(zkv_data$Estatura.ao.nascer[zkv_data$Group == "1"])
shapiro.test(zkv_data$Estatura.ao.nascer[zkv_data$Group == "0"])
assoc_heigh <- wilcox.test(Estatura.ao.nascer ~ Group, data = zkv_data)
print(assoc_heigh)

summary_heigh <- aggregate(Estatura.ao.nascer ~ Group, data = zkv_data,
                           FUN = function(x) c(median=median(x, na.rm=TRUE),
                                               Q1 = quantile(x, 0.25, na.rm=TRUE),                                 
                                               Q3 = quantile(x, 0.75, na.rm=TRUE),
                                               IQR=IQR(x, na.rm=TRUE),
                                               min=min(x, na.rm=TRUE), max=max(x, na.rm=TRUE)))
print(summary_heigh)

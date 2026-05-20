# Set working directory
setwd("~/Documents/Estatistica_mestrado")

# Load packages
library(openxlsx)
library(dplyr)
library(HardyWeinberg)
library(stringr)
library(SNPassoc)
library(haplo.stats)
library(textshape)
library(genetics)
library(snpStats)
library(LDheatmap)
library(viridis)
library(grid)
library(corrplot)
library(pwr)

# Load files
zkv_genotypes <- read.xlsx("Banco ZIKV - IL1B.xlsx")
zkv_metadata <- read.xlsx("metadados_zikv.xlsx")

# Prepare Data
zkv_data <- merge(zkv_metadata, zkv_genotypes, by = "Sample")

zkv_data <- zkv_data %>%
  mutate(
    Skin_color = case_when(
    `Raça/cor` == 1 ~ "White",
    `Raça/cor` == 2 ~ "Black",
    `Raça/cor` == 4 ~ "Brown"),
    
    Skin_color_recoded = case_when(
      `Raça/cor` == 1 ~ "White",
      `Raça/cor` %in% c(2, 4) ~ "Black/Brown"))


#___ Hardy-Weinberg Equilibrium ___#
run_HWE <- function(dataset, group, snp) {
  genotypes <- table(dataset[dataset$Group == group, snp])
  HWEex <- HWExact(genotypes, alternative = "two.sided", pvaluetype = "selome", eps=1e-10, x.linked = FALSE, verbose = TRUE)
}

HWEex_rs306_case <- run_HWE(zkv_data, 1, "rs4848306")
HWEex_rs306_ctrl <- run_HWE(zkv_data, 0, "rs4848306")
HWEex_rs623_case <- run_HWE(zkv_data, 1, "rs1143623")
HWEex_rs623_ctrl <- run_HWE(zkv_data, 0, "rs1143623")
HWEex_rs944_case <- run_HWE(zkv_data, 1, "rs16944")
HWEex_rs944_ctrl <- run_HWE(zkv_data, 0, "rs16944")
HWEex_rs627_case <- run_HWE(zkv_data, 1, "rs1143627")
HWEex_rs627_ctrl <- run_HWE(zkv_data, 0, "rs1143627")


#___ Allelic Frequency Comparison ___#
## Between groups
run_assoc_allele <- function(dataset, case, control, snp, allele1, allele2) {
  allele1_case <- sum(str_count(dataset[dataset$Group == case, snp], allele1))
  allele1_control <- sum(str_count(dataset[dataset$Group == control, snp], allele1))
  allele2_case <- sum(str_count(dataset[dataset$Group == case, snp], allele2))
  allele2_control <- sum(str_count(dataset[dataset$Group == control, snp], allele2))
  
  alleles <- matrix(c(allele1_case, allele2_case, allele1_control, allele2_control), nrow = 2, byrow = TRUE,
                    dimnames = list(Group = c("Case", "Control"), Allele = c(allele1, allele2)))
  
  assoc_alelos_test <- chisq.test(alleles)
  print(alleles)
  print(assoc_alelos_test)
}

assoc_alelos_rs306 <- run_assoc_allele(zkv_data, 1, 0, "rs4848306", "A", "G")
assoc_alelos_rs623 <- run_assoc_allele(zkv_data, 1, 0, "rs1143623", "C", "G")
assoc_alelos_rs944 <- run_assoc_allele(zkv_data, 1, 0, "rs16944", "A", "G")
assoc_alelos_rs627 <- run_assoc_allele(zkv_data, 1, 0, "rs1143627", "A", "G")

## vs. databases (ABraOM and DNABR)
compare_to_database <- function(dataset, snp, allele1, allele2, freq_db_allele1, freq_db_allele2) {
  allele1_case <- sum(str_count(dataset[dataset$Group == 1, snp], allele1))
  allele2_case <- sum(str_count(dataset[dataset$Group == 1, snp], allele2))
  allele1_control <- sum(str_count(dataset[dataset$Group == 0, snp], allele1))
  allele2_control <- sum(str_count(dataset[dataset$Group == 0, snp], allele2))
  
  print(paste0("Case: ", allele1, "=", allele1_case, "; ", allele2, "=", allele2_case))
  case_test <- chisq.test(c(allele1_case, allele2_case), p = c(freq_db_allele1, freq_db_allele2))
  print(case_test)
  
  print(paste0("Control: ", allele1, "=", allele1_control, "; ", allele2, "=", allele2_control))
  control_test <- chisq.test(c(allele1_control, allele2_control), p = c(freq_db_allele1, freq_db_allele2))
  print(control_test)
}

rs306_ABraOM <- compare_to_database(zkv_data, "rs4848306", "G", "A", 0.5918, 0.4082)
rs623_ABraOM <- compare_to_database(zkv_data, "rs1143623", "C", "G", 0.7447, 0.2553)
rs944_ABraOM <- compare_to_database(zkv_data, "rs16944", "A", "G", 0.4108, 0.5892)
rs627_ABraOM <- compare_to_database(zkv_data, "rs1143627", "A", "G", 0.5734, 0.4266)

rs306_DNABr <- compare_to_database(zkv_data, "rs4848306", "G", "A", 0.6375, 0.3625)
rs623_DNABr <- compare_to_database(zkv_data, "rs1143623", "C", "G", 0.7224, 0.2776)
rs944_DNABr <- compare_to_database(zkv_data, "rs16944", "A", "G", 0.4690, 0.5310)
rs627_DNABr <- compare_to_database(zkv_data, "rs1143627", "A", "G", 0.5075, 0.4925)

## Between self-reported skin color
run_assoc_allele_skin_color <- function(dataset, color1, color2, color3, snp, allele1, allele2) {
  allele1_color1 <- sum(str_count(dataset[dataset$Skin_color == color1, snp], allele1))
  allele1_color2 <- sum(str_count(dataset[dataset$Skin_color == color2, snp], allele1))
  allele1_color3 <- sum(str_count(dataset[dataset$Skin_color == color3, snp], allele1))
  allele2_color1 <- sum(str_count(dataset[dataset$Skin_color == color1, snp], allele2))
  allele2_color2 <- sum(str_count(dataset[dataset$Skin_color == color2, snp], allele2))
  allele2_color3 <- sum(str_count(dataset[dataset$Skin_color == color3, snp], allele2))
  
  alleles <- matrix(c(allele1_color1, allele1_color2, allele1_color3,
                      allele2_color1, allele2_color2, allele2_color3), nrow = 2, byrow = TRUE,
                    dimnames = list(Allele = c(allele1, allele2),
                                    Skin_color = c(color1, color2, color3)))
  
  assoc_alelos_test <- chisq.test(alleles)
  print(alleles)
  print(assoc_alelos_test)
}

assoc_skincolor_rs306 <- run_assoc_allele_skin_color(zkv_data, "White", "Black", "Brown", "rs4848306", "A", "G")
assoc_skincolor_rs623 <- run_assoc_allele_skin_color(zkv_data, "White", "Black", "Brown", "rs1143623", "C", "G")
assoc_skincolor_rs944 <- run_assoc_allele_skin_color(zkv_data, "White", "Black", "Brown", "rs16944", "A", "G")
assoc_skincolor_rs627 <- run_assoc_allele_skin_color(zkv_data, "White", "Black", "Brown", "rs1143627", "A", "G")

#___ Genotypic Association Analysis ___#
snpdata <- setupSNP(data = zkv_data, 
                    colSNPs = 20:23,  
                    sep = "",
                    info = zkv_data)

snpdata1 <- snpdata
snpdata1$rs16944  <- reorder(snpdata1$rs16944,  "minor")
snpdata1$rs1143627  <- reorder(snpdata1$rs1143627,  "minor")

association(Group ~ rs4848306, data = snpdata, model = "log-additive")
association(Group ~ rs1143623, data = snpdata, model = "log-additive")
association(Group ~ rs16944, data = snpdata1, model = "log-additive")
association(Group ~ rs1143627, data = snpdata1, model = "log-additive")

#___ Linkage Desequilibrium Analysis ___#
## Genetic R package
zkv_gt_snps <- zkv_genotypes %>%
  dplyr::select(-Sample)

zkv_gt_genetics <- zkv_gt_snps %>%
  mutate(rs1143627 = as.genotype(rs1143627, sep = ""),
         rs16944 = as.genotype(rs16944, sep = ""),
         rs1143623 = as.genotype(rs1143623, sep = ""),
         rs4848306 = as.genotype(rs4848306, sep = ""))

ld_genetics <- LD(zkv_gt_genetics)
print(ld_genetics)

## snpStats R package
zkv_gt_snstat <-  zkv_gt_snps %>%
  mutate(
    rs1143627 = case_when(
      rs1143627 == "GG" ~ 2,
      rs1143627 == "AG" ~ 1,
      rs1143627 == "AA" ~ 0),
    
    rs16944 = case_when(
      rs16944 == "AA" ~ 0,
      rs16944 == "AG" ~ 1,
      rs16944 == "GG" ~ 2),
    
    rs1143623 = case_when(
      rs1143623 == "CC" ~ 2,
      rs1143623 == "CG" ~ 1,
      rs1143623 == "GG" ~ 0),
    
    rs4848306 = case_when(
      rs4848306 == "GG" ~ 0,
      rs4848306 == "AG" ~ 1,
      rs4848306 == "AA" ~ 2))

snp_matrix <- as(as.matrix(zkv_gt_snstat), "SnpMatrix")
ld_matrix <- ld(snp_matrix, stats = c("D.prime", "R.squared"), depth = 3)
print(ld_matrix)

tiff("LD_r2.tiff", width = 5, height = 5, units = "in", res = 300)
LDheatmap(snp_matrix, 
          genetic.distances = c(-31, -511, -1464, -3737),
          LDmeasure = "r",
          SNP.name = c("rs1143627", "rs16944", "rs1143623", "rs4848306"),
          distances = "physical",
          add.map = TRUE,
          color = magma(20),
          title = NULL,
          geneMapLabelY = 0.2,
          name="ldheatmap")
grid.edit(gPath("ldheatmap", "geneMap","SNPnames"), gp = gpar(cex=1, col="black"))
dev.off()

#___ Haplotype Association Analysis ___#
snpsH <- c("rs4848306", "rs1143623",  "rs16944", "rs1143627")
genoH <- make.geno(snpdata, snpsH)
print(genoH)

haplo_cc_results <- haplo.cc(
  y = snpdata$Group,
  geno = genoH,
  x.adj = NULL,
  locus.label = snpsH,
  ci.prob = 0.95,
  simulate = FALSE)

print(haplo_cc_results)

p_vals <- haplo_cc_results[["cc.df"]][["p-val"]]
p_vals_adjusted <- p.adjust(p_vals, method = "bonferroni")


#___ Logistic Regression Analysis ___# 
hap_assignments <- haplo.em(genoH, locus.label = snpsH) %>% 
  summary() %>% 
  as.data.frame() %>% 
  dplyr::filter(posterior > 0.90) %>%
  dplyr::mutate(ID = snpdata$Sample[subj.id],
                Group = snpdata$Group[subj.id],
                Carrier = (hap1 == 4 | hap2 == 4))

zkv_data <- zkv_data %>%
  left_join(hap_assignments %>% dplyr::select(ID, hap1, hap2, posterior, Carrier), by = c("Sample" = "ID"))

## Unadjusted (haplotype only)
model_unadj <- glm(Group ~ Carrier, data = zkv_data, family = binomial)

summary(model_unadj)
or_unadj <- exp(coef(model_unadj))
ci_unadj <- exp(confint(model_unadj))

## Adjusted for skin color
model_skin <- glm(Group ~ Carrier + Skin_color_recoded, data = zkv_data, family = binomial)

summary(model_skin)
or_skin <- exp(coef(model_skin))
ci_skin <- exp(confint(model_skin))

## Adjusted for trimester
zkv_complete_trim <- zkv_data %>% filter(!is.na(Trim..Sint..de.Zika)) # complete cases
print(table(zkv_complete_trim$Group))

model_trim <- glm(Group ~ Carrier + factor(Trim..Sint..de.Zika), data = zkv_complete_trim, family = binomial)

summary(model_trim)
or_trim <- exp(coef(model_trim))
ci_trim <- exp(confint(model_trim))

## Full model (trimester + skin color)
model_full <- glm(Group ~ Carrier + factor(Trim..Sint..de.Zika) + Skin_color_recoded, data = zkv_complete_trim, family = binomial)

summary(model_full)
or_full <- exp(coef(model_full))
ci_full <- exp(confint(model_full))

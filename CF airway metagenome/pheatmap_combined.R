library(readr)
library(readxl)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(compositions)
library(ggpubr)

RPMM_commensals <- read_delim("CF_RPMM_commensals.csv", 
                              delim = ";", escape_double = FALSE, trim_ws = TRUE)
RPMM_commensals <- data.frame(RPMM_commensals)
rownames(RPMM_commensals) <- RPMM_commensals$species
RPMM_commensals$species <- NULL



RPMM_HI <- read_delim("CF_RPMM_HI.csv", 
                      delim = ";", escape_double = FALSE, trim_ws = TRUE)
RPMM_HI <- data.frame(RPMM_HI)
rownames(RPMM_HI) <- RPMM_HI$species
RPMM_HI$species <- NULL


RPMM_PA_SA <- read_delim("CF_RPMM_PA_SA.csv", 
                         delim = ";", escape_double = FALSE, trim_ws = TRUE)
RPMM_PA_SA <- data.frame(RPMM_PA_SA)
rownames(RPMM_PA_SA) <- RPMM_PA_SA$species
RPMM_PA_SA$species <- NULL


RPMM_PA <- read_delim("CF_RPMM_PA.csv", 
                      delim = ";", escape_double = FALSE, trim_ws = TRUE)
RPMM_PA <- data.frame(RPMM_PA)
rownames(RPMM_PA) <- RPMM_PA$species
RPMM_PA$species <- NULL

RPMM_SA <- read_delim("CF_RPMM_SA.csv", 
                      delim = ";", escape_double = FALSE, trim_ws = TRUE)
RPMM_SA <- data.frame(RPMM_SA)
rownames(RPMM_SA) <- RPMM_SA$species
RPMM_SA$species <- NULL


RPMM_PA_SM <- read_delim("CF_RPMM_PA_SM.csv", 
                         delim = ";", escape_double = FALSE, trim_ws = TRUE)
RPMM_PA_SM <- data.frame(RPMM_PA_SM)
rownames(RPMM_PA_SM) <- RPMM_PA_SM$species
RPMM_PA_SM$species <- NULL


RPMM_SM_SA <- read_delim("CF_RPMM_SM_SA.csv", 
                         delim = ";", escape_double = FALSE, trim_ws = TRUE)
RPMM_SM_SA <- data.frame(RPMM_SM_SA)
rownames(RPMM_SM_SA) <- RPMM_SM_SA$species
RPMM_SM_SA$species <- NULL


RPMM_SM <- read_delim("CF_RPMM_SM.csv", 
                      delim = ";", escape_double = FALSE, trim_ws = TRUE)
RPMM_SM <- data.frame(RPMM_SM)
rownames(RPMM_SM) <- RPMM_SM$species
RPMM_SM$species <- NULL


RPMM_healthy <- read_excel("RPMM_healthy.xlsx")
RPMM_healthy <- data.frame(RPMM_healthy)
rownames(RPMM_healthy) <- RPMM_healthy$species
RPMM_healthy$species <- NULL


#commensals
RPMM_commensals_00 <- merge(RPMM_healthy, RPMM_commensals,
                            by = "row.names", all.x = TRUE, all.y = TRUE)
rownames(RPMM_commensals_00) <- RPMM_commensals_00$Row.names
RPMM_commensals_00$Row.names <- NULL

RPMM_commensals_00 <- RPMM_commensals_00[order(rowSums(RPMM_commensals_00), decreasing = TRUE),]

RPMM_commensals_01 <- RPMM_commensals_00[c(1:50),]

RPMM_commensals_01 <- RPMM_commensals_01 %>% 
  mutate_if(is.numeric, ~round(., 0))
RPMM_commensals_01[RPMM_commensals_01 == 0] <- NA


pheatmap(log10(RPMM_commensals_01+1),
         na_col = "grey40",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = 1,
         color=colorRampPalette(c("black", "blue", "red","yellow"))(10000))

#HI
RPMM_HI_00 <- merge(RPMM_healthy, RPMM_HI,
                            by = "row.names", all.x = TRUE, all.y = TRUE)
rownames(RPMM_HI_00) <- RPMM_HI_00$Row.names
RPMM_HI_00$Row.names <- NULL

RPMM_HI_00 <- RPMM_HI_00[order(rowSums(RPMM_HI_00), decreasing = TRUE),]

RPMM_HI_01 <- RPMM_HI_00[c(1:50),]

RPMM_HI_01 <- RPMM_HI_01 %>% 
  mutate_if(is.numeric, ~round(., 0))
RPMM_HI_01[RPMM_HI_01 == 0] <- NA


pheatmap(log10(RPMM_HI_01+1),
         na_col = "grey40",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = 1,
         color=colorRampPalette(c("black", "blue", "red","yellow"))(10000))


#SM
RPMM_SM_00 <- merge(RPMM_healthy, RPMM_SM,
                    by = "row.names", all.x = TRUE, all.y = TRUE)
rownames(RPMM_SM_00) <- RPMM_SM_00$Row.names
RPMM_SM_00$Row.names <- NULL
RPMM_SM_00[is.na(RPMM_SM_00)] <- 0

RPMM_SM_001 <- merge(RPMM_SM_00, RPMM_PA_SM,
                     by = "row.names", all.x = TRUE, all.y = TRUE)
rownames(RPMM_SM_001) <- RPMM_SM_001$Row.names
RPMM_SM_001$Row.names <- NULL
RPMM_SM_001[is.na(RPMM_SM_001)] <- 0

RPMM_SM_002 <- merge(RPMM_SM_001, RPMM_SM_SA,
                     by = "row.names", all.x = TRUE, all.y = TRUE)
rownames(RPMM_SM_002) <- RPMM_SM_002$Row.names
RPMM_SM_002$Row.names <- NULL
RPMM_SM_002[is.na(RPMM_SM_002)] <- 0


RPMM_SM_002 <- RPMM_SM_002[order(rowSums(RPMM_SM_002), decreasing = TRUE),]

RPMM_SM_01 <- RPMM_SM_002[c(1:50),]

RPMM_SM_01 <- RPMM_SM_01 %>% 
  mutate_if(is.numeric, ~round(., 0))
RPMM_SM_01[RPMM_SM_01 == 0] <- NA


pheatmap(log10(RPMM_SM_01+1),
         na_col = "grey40",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = 1,
         color=colorRampPalette(c("black", "blue", "red","yellow"))(10000))

#combine commensals, HI, SM to one big heatmap

RPMM_comHISM <- merge(RPMM_commensals_00, RPMM_HI_00,
                      by = "row.names", all.x = TRUE, all.y = TRUE)
rownames(RPMM_comHISM) <- RPMM_comHISM$Row.names
RPMM_comHISM$Row.names <- NULL
RPMM_comHISM[is.na(RPMM_comHISM)] <- 0

RPMM_comHISM <- merge(RPMM_comHISM, RPMM_SM_002,
                      by = "row.names", all.x = TRUE, all.y = TRUE)
rownames(RPMM_comHISM) <- RPMM_comHISM$Row.names
RPMM_comHISM$Row.names <- NULL
RPMM_comHISM[is.na(RPMM_comHISM)] <- 0

RPMM_comHISM <- RPMM_comHISM[order(rowSums(RPMM_comHISM), decreasing = TRUE),]

RPMM_comHISM_01 <- RPMM_comHISM[c(1:50),]

RPMM_comHISM_01 <- RPMM_comHISM_01 %>% 
  mutate_if(is.numeric, ~round(., 0))
RPMM_comHISM_01[RPMM_comHISM_01 == 0] <- NA


comHISM <- pheatmap(log10(RPMM_comHISM_01+1),
         na_col = "grey40",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = c(1,15,16,20,21),
         color=colorRampPalette(c("black", "blue", "red","yellow"))(10000))

pdf("comHISM.pdf", width=10, height=8)

comHISM

dev.off()


#PA
RPMM_PA_00 <- merge(RPMM_healthy, RPMM_PA,
                            by = "row.names", all.x = TRUE, all.y = TRUE)
rownames(RPMM_PA_00) <- RPMM_PA_00$Row.names
RPMM_PA_00$Row.names <- NULL
RPMM_PA_00[is.na(RPMM_PA_00)] <- 0

RPMM_PA_00 <- RPMM_PA_00[order(rowSums(RPMM_PA_00), decreasing = TRUE),]

RPMM_PA_01 <- RPMM_PA_00[c(1:50),]

RPMM_PA_01 <- RPMM_PA_01 %>% 
  mutate_if(is.numeric, ~round(., 0))
RPMM_PA_01[RPMM_PA_01 == 0] <- NA

pheatmap(log10(RPMM_PA_01+1),
         na_col = "grey40",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = c(1),
         color=colorRampPalette(c("black", "blue", "red","yellow"))(10000))

#SA
RPMM_SA_00 <- merge(RPMM_healthy, RPMM_SA,
                    by = "row.names", all.x = TRUE, all.y = TRUE)
rownames(RPMM_SA_00) <- RPMM_SA_00$Row.names
RPMM_SA_00$Row.names <- NULL
RPMM_SA_00[is.na(RPMM_SA_00)] <- 0

RPMM_SA_00 <- RPMM_SA_00[order(rowSums(RPMM_SA_00), decreasing = TRUE),]

RPMM_SA_01 <- RPMM_SA_00[c(1:50),]

RPMM_SA_01 <- RPMM_SA_01 %>% 
  mutate_if(is.numeric, ~round(., 0))
RPMM_SA_01[RPMM_SA_01 == 0] <- NA

pheatmap(log10(RPMM_SA_01+1),
         na_col = "grey40",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = c(1),
         color=colorRampPalette(c("black", "blue", "red","yellow"))(10000))

#PA_SA
RPMM_PA_SA_00 <- merge(RPMM_healthy, RPMM_PA_SA,
                    by = "row.names", all.x = TRUE, all.y = TRUE)
rownames(RPMM_PA_SA_00) <- RPMM_PA_SA_00$Row.names
RPMM_PA_SA_00$Row.names <- NULL
RPMM_PA_SA_00[is.na(RPMM_PA_SA_00)] <- 0

RPMM_PA_SA_00 <- RPMM_PA_SA_00[order(rowSums(RPMM_PA_SA_00), decreasing = TRUE),]

RPMM_PA_SA_01 <- RPMM_PA_SA_00[c(1:50),]

RPMM_PA_SA_01 <- RPMM_PA_SA_01 %>% 
  mutate_if(is.numeric, ~round(., 0))
RPMM_PA_SA_01[RPMM_PA_SA_01 == 0] <- NA

pheatmap(log10(RPMM_PA_SA_01+1),
         na_col = "grey40",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = c(1),
         color=colorRampPalette(c("black", "blue", "red","yellow"))(10000))


#combine PA, SA, PA_SA to one big heatmap

RPMM_PASA <- merge(RPMM_PA_00, RPMM_SA_00,
                      by = "row.names", all.x = TRUE, all.y = TRUE)
rownames(RPMM_PASA) <- RPMM_PASA$Row.names
RPMM_PASA$Row.names <- NULL
RPMM_PASA[is.na(RPMM_PASA)] <- 0

RPMM_PASA <- merge(RPMM_PASA, RPMM_PA_SA_00,
                      by = "row.names", all.x = TRUE, all.y = TRUE)
rownames(RPMM_PASA) <- RPMM_PASA$Row.names
RPMM_PASA$Row.names <- NULL
RPMM_PASA[is.na(RPMM_PASA)] <- 0

RPMM_PASA <- RPMM_PASA[order(rowSums(RPMM_PASA), decreasing = TRUE),]

RPMM_PASA_01 <- RPMM_PASA[c(1:50),]

RPMM_PASA_01 <- RPMM_PASA_01 %>% 
  mutate_if(is.numeric, ~round(., 0))
RPMM_PASA_01[RPMM_PASA_01 == 0] <- NA


PASA <- pheatmap(log10(RPMM_PASA_01+1),
         na_col = "grey40",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = c(1,16,17,32,33),
         color=colorRampPalette(c("black", "blue", "red","yellow"))(10000))

pdf("PASA.pdf", width=10, height=8)

PASA

dev.off()


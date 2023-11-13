setwd("/mnt/sfb900nfs/groups/tuemmler/ilona/nina_BE/R_paper")

library(readxl)
library(ggpubr)
library(RColorBrewer)

final_table <- read_excel("neu_MLST_hexCode_merged.xlsx")
final_table <- data.frame(final_table)

final_table <- subset(final_table, MLST != "-")

counts_MHH <- table(final_table$MLST)
counts_MHH <- data.frame(counts_MHH)

colnames(counts_MHH) <- c("MLST", "count")
rownames(counts_MHH) <- counts_MHH$MLST

counts_MHH$MLST <- as.numeric(as.character(counts_MHH$MLST))
counts_MHH <- counts_MHH[order(counts_MHH$MLST),]
counts_MHH$MLST <- as.character(counts_MHH$MLST)

counts_MHH$group <- with(counts_MHH,
                         ifelse(counts_MHH$MLST == 395, "#FDB462",
                         ifelse(counts_MHH$MLST == 253, "#8DD3C7",
                         ifelse(counts_MHH$MLST == 27, "#BEBADA",
                         ifelse(counts_MHH$MLST == 17, "#FB8072",
                         #ifelse(counts_MHH$MLST == 274, "#80B1D3", 
                         "black")))))
#)

display.brewer.pal(6, name = "Set3")
brewer.pal(6, name = "Set3")


ggdotchart(counts_MHH, x = "MLST", y = "count",
           color = "group",
           palette = c(#"#80B1D3", 
             "#8DD3C7", "#BEBADA", "#FB8072", "#FDB462", "grey40"),
           sorting = "none",
           add = "segments",
           dot.size = 6,
           label = round(counts_MHH$count),
           font.label = list(color = "white", size = 9, vjust = 0.5),
           add.params = list(color = "group", size = 1),
           x.text.col = TRUE,
           ylab = "Count",
           xlab = "Sequence types",
           legend = "none",
           ggtheme = theme_pubr()) +
  scale_y_continuous(limits = c(0,10),
                     breaks = c(0,5,10)) +
  scale_x_discrete(expand = c(0.02,0)) +
  theme(axis.text.x = element_text(size = 10, face = "bold"))


#human database
human <- read_excel("pseudomonas_database_filtered_R.xlsx", 
                    sheet = "filtered_human")
human <- data.frame(human)

counts_human <- table(human$MLST)
counts_human <- data.frame(counts_human)

colnames(counts_human) <- c("MLST", "count")
rownames(counts_human) <- counts_human$MLST

counts_human$MLST <- as.numeric(as.character(counts_human$MLST))
counts_human <- counts_human[order(counts_human$MLST),]
counts_human$MLST <- as.character(counts_human$MLST)

ggdotchart(counts_human, x = "MLST", y = "count",
           sorting = "descending",
           add = "segments",
           dot.size = 4,
           label = round(counts_human$count),
           font.label = list(color = "white", size = 7, vjust = 0.5),
           ggtheme = theme_pubr()) +
  scale_y_continuous(limits = c(0,45),
                     breaks = c(0,20,40)) +
  scale_x_discrete(expand = c(0.005,0)) +
  theme(axis.text.x = element_text(size=8))

##
counts_human <- counts_human[order(counts_human$count, decreasing = TRUE),]
counts_human$relAbund <- round(counts_human$count / sum(counts_human$count) * 100, 1)
counts_human_01 <- counts_human[c(1:19),]
counts_human_01$ST <- "ST"
counts_human_01$ST <- paste(counts_human_01$ST, counts_human_01$MLST, sep ="")


ggplot(counts_human_01, aes(x = reorder(ST, -count), y = relAbund, fill = ST)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(breaks = c ("ST235", "ST111", "ST17", "ST179", "ST308", "ST146", "ST244",
                                "ST27", "ST155", "ST348", "ST274", "ST277", "ST253", "ST395",
                                "ST845", "ST245", "ST170", "ST242", "ST649"),
                    values = c("grey40", "grey40", "#FB8072", "grey40", "grey40", "grey40", "grey40",
                               "#BEBADA", "grey40", "grey40", "#80B1D3", "grey40", "#8DD3C7", "#FDB462", 
                               "grey40", "grey40", "grey40", "grey40", "grey40")) +
  scale_y_continuous(limits = c(0,6.5),
                     breaks = c(0,2,4,6),
                     expand = c(0,0.1)) +
  ylab("Relative abundance (%)") +
  xlab("Human infection") +
  guides(fill = FALSE) +
  theme_pubr()

#human chronic database
chronic <- read_excel("pseudomonas_database_filtered_R.xlsx", 
                    sheet = "filtered_chronic")
chronic <- data.frame(chronic)

counts_chronic <- table(chronic$MLST)
counts_chronic <- data.frame(counts_chronic)

colnames(counts_chronic) <- c("MLST", "count")
rownames(counts_chronic) <- counts_chronic$MLST

counts_chronic$MLST <- as.numeric(as.character(counts_chronic$MLST))
counts_chronic <- counts_chronic[order(counts_chronic$MLST),]
counts_chronic$MLST <- as.character(counts_chronic$MLST)

ggdotchart(counts_chronic, x = "MLST", y = "count",
           sorting = "descending",
           add = "segments",
           dot.size = 4,
           label = round(counts_chronic$count),
           font.label = list(color = "white", size = 7, vjust = 0.5),
           ggtheme = theme_pubr()) +
  scale_y_continuous(limits = c(0,35),
                     breaks = c(0,10,30)) +
  scale_x_discrete(expand = c(0.005,0)) +
  theme(axis.text.x = element_text(size=8))

##
counts_chronic <- counts_chronic[order(counts_chronic$count, decreasing = TRUE),]
counts_chronic$relAbund <- round(counts_chronic$count / sum(counts_chronic$count) * 100, 1)
counts_chronic_01 <- counts_chronic[c(1:19),]
counts_chronic_01$ST <- "ST"
counts_chronic_01$ST <- paste(counts_chronic_01$ST, counts_chronic_01$MLST, sep ="")


ggplot(counts_chronic_01, aes(x = reorder(ST, -count), y = relAbund, fill = ST)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(breaks = c ("ST17", "ST146", "ST179", "ST27", "ST155", "ST111", "ST274",
                                "ST244", "ST170", "ST395", "ST649", "ST242", "ST253", "ST348",
                                "ST132", "ST471", "ST775", "ST782", "ST845"),
                    values = c("#FB8072", "grey40", "grey40", "#BEBADA", "grey40", "grey40", "#80B1D3",
                               "grey40", "grey40", "#FDB462", "grey40", "grey40", "#8DD3C7", "grey40",
                               "grey40", "grey40", "grey40", "grey40", "grey40")) +
  scale_y_continuous(limits = c(0,6.5),
                     breaks = c(0,2,4,6),
                     expand = c(0,0.1)) +
  ylab("Relative abundance (%)") +
  xlab("Chronic human infection") +
  guides(fill = FALSE) +
  theme_pubr() +
  theme(axis.text.x = element_text(size = 10))

#environment
environment <- read_excel("pseudomonas_database_filtered_R.xlsx", 
                      sheet = "filtered_environment")
environment <- data.frame(environment)

counts_environment <- table(environment$MLST)
counts_environment <- data.frame(counts_environment)

colnames(counts_environment) <- c("MLST", "count")
rownames(counts_environment) <- counts_environment$MLST

counts_environment$MLST <- as.numeric(as.character(counts_environment$MLST))
counts_environment <- counts_environment[order(counts_environment$MLST),]
counts_environment$MLST <- as.character(counts_environment$MLST)

#hospital environment
hospital <- read_excel("pseudomonas_database_filtered_R.xlsx", 
                          sheet = "filtered_hospital")
hospital <- data.frame(hospital)

counts_hospital <- table(hospital$MLST)
counts_hospital <- data.frame(counts_hospital)

colnames(counts_hospital) <- c("MLST", "count")
rownames(counts_hospital) <- counts_hospital$MLST

counts_hospital$MLST <- as.numeric(as.character(counts_hospital$MLST))
counts_hospital <- counts_hospital[order(counts_hospital$MLST),]
counts_hospital$MLST <- as.character(counts_hospital$MLST)


#export tables with counts

MHH <- counts_MHH
colnames(MHH) <- c("MLST", "counts_MHH")

human <- counts_human
colnames(human) <- c("MLST", "counts_human")

chronic <- counts_chronic
colnames(chronic) <- c("MLST", "counts_chronic")

environment <- counts_environment
colnames(environment) <- c("MLST", "counts_environment")

hospital <- counts_hospital
colnames(hospital) <- c("MLST", "counts_hospital")

all <- merge(MHH, human, by = 0, all = TRUE)

rownames(all) <- all$Row.names
all$Row.names <- NULL
all$MLST.x <- NULL
all$MLST.y <- NULL

all <- merge(all, chronic, by = 0, all = TRUE)

rownames(all) <- all$Row.names
all$Row.names <- NULL
all$MLST <- NULL

all <- merge(all, environment, by = 0, all = TRUE)

rownames(all) <- all$Row.names
all$Row.names <- NULL
all$MLST <- NULL

all <- merge(all, hospital, by = 0, all = TRUE)

rownames(all) <- all$Row.names
all$Row.names <- NULL
all$MLST <- NULL

all$MLST <- rownames(all)
all$MLST <- as.numeric(all$MLST)
all <- all[order(all$MLST),]

all[is.na(all)] <- 0

write.table(all, file = "MLST_MHH_human_chronic.csv", sep = ";",
            row.names = TRUE)

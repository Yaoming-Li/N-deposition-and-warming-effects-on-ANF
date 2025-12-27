library(openxlsx)
library(agricolae)
library(ggplot2)

############variance analysis###########
Env = read.xlsx('Env.xlsx',1,rowNames = TRUE)
D<-aov(Env$ANFR~Env$Treatment,data = Env)
anova(D)

## Duncan's Multiple Comparisons
P <-HSD.test(D,"Env$Treatment")

## Calculate the mean and standard deviation of each group
mar<-P$groups
rownamemar<-row.names(mar)
newmar<-data.frame(rownamemar,mar[1],mar$groups)
sort<-newmar[order(newmar$rownamemar),]
Treatment<-row.names(P$means)
mean<-P$means[,1]
sd<-P$means[,2]
marker<-sort$mar.groups
plotdata<-data.frame(Treatment,mean,sd,marker)

#boxplot
idcol<-(c("gray","#A7D3D4","#009B9E", "#E4C1D9","#C75DAB"))
p <- ggplot(Env, aes(x = Treatment, y = ANFR,fill=Treatment)) + 
  stat_boxplot(geom = "errorbar", width = 0.15, size = 0.5) +
  geom_boxplot(width = 0.4, alpha = 1) +
  stat_summary(fun = "mean", geom = "point", shape = 8, size = 3) +
 
  theme(strip.text.x = element_text(size = 20, angle = 0)) +  
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  labs(x = '', y = 'ANFR') +
  theme_bw() +  
  theme(panel.grid = element_blank(),  
        axis.ticks.length.x = unit(0.1, "cm"),  
        axis.ticks.length.y = unit(0.1, "cm"), 
        axis.line = element_line(color = "black"),  
        axis.text = element_text(size = 12, color = "black")) +  
  scale_fill_manual(values = idcol) + 
   geom_text(data = plotdata, aes(x = Treatment, y = max(Env$ANFR)+0.5, label = marker), 
            size = 5, color = "black")


##################Calculate the effect size##############
library(tidyverse)
library(metafor)
#env <- as.data.frame(scale(Env[,-c(1,2)]))
#env<-cbind(env+2,Env[,c(1,2)])

#  Calculate the mean and standard deviation of each group
summary_data <- Env %>%  filter(Treatment != "CK") %>% 
  group_by(Treatment) %>%  summarize(
    TMean = mean(TC),
    TSD = sd(TC),
    TN = n(), 
    CMean = mean(Env$TC[Env$Treatment == "CK"]),      
    CSD = sd(Env$TC[Env$Treatment == "CK"]),  
    C_NP = n()   
  )

d1 <- summary_data
d2 <- escalc(measure = "ROM", data = d1, m1i = TMean, sd1i = TSD, n1i = TN, m2i = CMean, sd2i = CSD, n2i = C_NP)
r1 <- rma(yi, vi, data = d2, method = "REML")
kk1 <- summary(r1)

subgroups <- unique(d2$Treatment) 
subgroup_results <- list()  
for (subgroup in subgroups) {
  subgroup_data <- subset(d2, Treatment == subgroup)    
  subgroup_result <- rma(yi, vi, data = subgroup_data, method = "REML")   
  subgroup_results[[subgroup]] <- subgroup_result}


for (subgroup in subgroups){ 
  print(subgroup_results[[subgroup]])
}

combined_data1 <- data.frame(estimate = kk1$b,se = kk1$se,pval = kk1$pval)
combined_subgroup_data <- list()
for (subgroup in subgroups) {  
  subgroup_result <- subgroup_results[[subgroup]]  
  combined_subgroup_data[[subgroup]] <- data.frame(
    estimate = subgroup_result$b,
    se = subgroup_result$se,
    pval = subgroup_result$pval  )}

combined_data_subgroup <- do.call(rbind, combined_subgroup_data)

combined_data_subgroup$Subgroup <- rownames(combined_data_subgroup)

data_df <- combined_data_subgroup
data_df$sig <- ""

for (i in 1:nrow(data_df)) {
  if (data_df$pval[i] < 0.001) {
    data_df$sig[i] <- "***"  
  } else if (data_df$pval[i] < 0.01) {    
    data_df$sig[i] <- "**"  
  } else if (data_df$pval[i] < 0.05) {   
    data_df$sig[i] <- "*"  }}
print(data_df)

library(ggthemes)

dataset <- as.data.frame(data_df)
dataset$index <- rownames(dataset)

dataset$index = factor(dataset$index, 
                       levels=c("N1", "N5", "W", "WN"))
col<-(c("#A7D3D4","#009B9E", "#E4C1D9",  "#C75DAB"))
q1 <- ggplot(dataset, aes(x = index, y = estimate)) +  
  geom_point(aes(color = index), size=4) +  
  geom_errorbar(aes(ymax = estimate + 1.96*se, ymin = estimate - 1.96*se, color = index), 
                width=0.2, size=0.5) +  
  geom_hline(aes(yintercept = 0), color="gray", linetype="dashed", size = 0.6) +  
  xlab('') +  ylab('TC Effect size') +  
  theme_few() +  
  geom_text(aes(x = index, y = estimate + 1.96*se + 0.05, label = sig), 
            size = 6, position = position_dodge(1)) +
  scale_y_continuous(breaks = scales::breaks_extended(n = 8), 
                                          expand = expansion(mult = c(0.15, 0.15))
  ) +
  theme(axis.text.x = element_text(size = 24, color = "black")) +  
  theme(axis.text.y = element_text(size = 24, color = "black")) +  
  theme(title = element_text(size = 24)) +  
  theme_bw() +  
  theme(axis.ticks.length = unit(-0.25, "cm")) +
  scale_color_manual(values = col) +  scale_fill_manual(values = col) + 
  theme(legend.position = "NA")

###Merge ASV and group information
library(ggpubr)
library(vegan)
otu = read.xlsx('Diazotroph.xlsx',1,rowNames = TRUE)
otu <- otu[rowSums(otu != 0) > 0, ]
OTU <-t(otu)
betadis <- merge(OTU, Env, by = "row.names")
rownames(betadis) <- betadis[,1]
betadis[,1] <- NULL

####NMDS analysis
nmds <- metaMDS(betadis[,-c(1031:1053)])
scrs <-scores(nmds, display = 'sites')

scrs <-cbind(as.data.frame(scrs),Group = betadis$Treatment)

cent <-aggregate(cbind(NMDS1, NMDS2) ~ Group, data = scrs, FUN = mean)
segs <-merge(scrs, setNames(cent, c('Group','oNMDS1','ONMDS2')),
             by ='Group',sort = FALSE)
library("ggalt")
idcol<-(c("gray","#A7D3D4","#009B9E", "#E4C1D9",  "#C75DAB"))
p<-ggplot(scrs, aes(x=NMDS1,y=NMDS2, colour = Group))+
    geom_point(size=2)+
  geom_polygon(alpha = 0.1) +
    annotate("text", x = -1.5, y = 0.1, 
           label = paste0("stress: ",format(nmds$stress, digits = 4)), hjust = 0) +
  scale_color_manual(values = idcol)+theme_bw() + 
  theme(panel.grid = element_blank(),  
               axis.line = element_line(color = "black"), 
        axis.text = element_text(size = 12, color = "black"))

#### Micro Permanova
adon <- adonis2(betadis[,-c(1031:1053)] ~ Treatment, data=betadis, permutations=999)

#################Calculate diversity index##############
library(picante) 
library(iCAMP)

shannon <- diversity(OTU,MARGIN = 1) 
simpson <- diversity(OTU,index = "simpson",MARGIN = 1)   
SR <- specnumber(OTU,MARGIN = 1)  
Pielou <- shannon/log(SR)    
Chao1 <- estimateR(OTU)[2, ]

tree <- read.tree('tree.nwk')

Faith_PD <- pd(OTU, tree, include.root = FALSE)
faith_pd <-Faith_PD$PD

### beta dispersal
library("ggpubr")
df <- t(decostand(OTU,"total")) 
dis <- vegdist(t(df), method = "bray", na.rm = T)
Treatment <- factor(Env$Treatment)
mod <- betadisper(dis, Treatment)
distan <- mod$distances
diversity <- data.frame(shannon,simpson,SR,Pielou,Chao1,faith_pd,distan)
diversity <- cbind(diversity, Treatment)


##########Test of the negative binomial distribution
###Obtain the p-value and logFC
library(DESeq2)

####The absolute abundance table is obtained by multiplying 
#the relative abundance table by the absolute abundance of the nifH gene.
micro_abab <- read.xlsx('Absolute_abundance.xlsx',1,rowNames = TRUE)
sample_info <-Env[,1,drop=FALSE]

#The grouping variable (such as "Treatment") needs to be of the "factor" type.
sample_info$Treatment <- factor(sample_info$Treatment)

# Adopt the official recommended standard of DESeq2 (rowSums >= 10)
otu_table <- micro_abab[rowSums(micro_abab) >= 10, ]

otu_table <- as.matrix(otu_table)
storage.mode(otu_table) <- "integer"  

dds <- DESeqDataSetFromMatrix(countData = otu_table,
                              colData = sample_info,
                              design = ~ Treatment) 
# Conduct a difference analysis
dds <- DESeq(dds)
keep <- rowSums(counts(dds) >= 10) >= 3 
dds <- dds[keep, ]

# Obtain the results of the difference analysis
results1 <- results(dds,contrast = c("Treatment", "N1", "CK"), alpha = 0.05)

results_N1 <- as.data.frame(results1)

Tax<- merge(results_N1, tax, by = "row.names")
rownames(Tax) <- Tax[,1]
Tax[,1] <- NULL

library(dplyr)
Tax <- Tax %>%
  mutate(change = case_when(
    padj < 0.05 & abs(log2FoldChange) > 1 ~ ifelse(log2FoldChange < 0, "depleted", "enriched"),
    TRUE ~ "nosig"))
write.csv(Tax, "differential ASV.csv")

###################differential ASV analysis################
# Classification based on p-value and logFC
library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)
library(RColorBrewer)
library(CMplot)

Tax <- read.xlsx('differential ASV.xlsx',1,rowNames = TRUE)

#Take the negative logarithm of the p-value
Tax <- Tax %>%mutate(padj = -log10(padj))
Tax$Genus[is.na(Tax$Genus)] <- "Others"

Genus_order <- c(" g__Immundisolibacter"," g__Mesorhizobium"," g__Burkholderia",
                 " g__Methylocystis"," g__Paenibacillus"," g__Pseudacidovorax",
                 " g__Skermanella", " g__Clostridium", " g__Bradyrhizobium",
                 " g__Desulfovibrio", " g__Nostoc"," g__Rhizobium","Others")
Tax$Genus <- factor(Tax$Genus, levels = Genus_order)

paired_colors <- brewer.pal(n = 8, name = "Dark2")
color_values <- colorRampPalette(paired_colors)(13)##扩展颜色数量
color_values <-c("#d5231d","#3777ac","#4ea64a","#8e4c99","#e88f18",
                 "#e47faf","#b698c5","#a05528","#58a6d6","#1f2d6f",
                 "#279772","#add387","#d9b71a")

###Manhattan Map##
p1 <- ggplot(Tax, aes(x = Genus, y = padj)) +
  geom_point(aes(color = Genus, shape = change, size = abs(log2FoldChange)), position = position_jitter(width = 0.5), alpha = 0.7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(-2, 15)) +
  scale_size_continuous(range = c(1, 3.5)) +
  scale_shape_manual(values = c(6,17,16)) +
  scale_color_manual(values = color_values) + #自定义color_values
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "gray") +
  labs(x = "Genus", y = "-log10(padj)", 
       title = "Differential ASV in NX vs CK") +
  theme_light() +
  theme(legend.position = "top", 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###Summary of Differences in ASV
data <- data.frame(
  Treatment = c("N1", "N5", "W", "WN"),
  enriched = c(14, 15, 13, 15),
  depleted = c(46, 11, 41, 55))
library(reshape2)
data_long <- melt(data, id.vars = "Treatment", variable.name = "Status", value.name = "Count")
col <- c("enriched" = "#d9b71a", "depleted" = "#b698c5")

# Draw a stacked bar chart
ggplot(data_long, aes(x = Treatment, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "stack") +  
  scale_fill_manual(values = col) +  
  labs(x = "Treatment", y = "Count", title = "Enriched and Depleted ASVs by Treatment") +  # 添加标题和轴标签
  theme_minimal() +  
  scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) + 
  labs( x = "Treatment", y = "Differential ASV") +
  theme_test() + theme(legend.position = "bottom")

#################Evolutionary tree##############
###removing the ASVs that did not show significant differences among the four treatments
library(ggnewscale)
library(ggplot2)
library(ggtreeExtra)
library(ggtree)
library(iCAMP)
library(picante)
logFC <- read.xlsx('differential ASV.xlsx',6,rowNames = TRUE)
tax = read.xlsx('Diazotroph.xlsx',3,rowNames = TRUE)
tree <- read.tree('tree.nwk')
spid.check=match.name(cn.list=list(otu=t(logFC)),tree.list=list(tree=tree))
tree=spid.check$tree
otu_tax = as.data.frame(merge(logFC,tax,by = "row.names",all = F))
colnames(otu_tax)[1] = c("ASV")
otu_tax[is.na(otu_tax) | otu_tax == " g__"] <- "Others"
groupInfo <- split(otu_tax$ASV, otu_tax$Genus) 
Tree <- groupOTU(tree, groupInfo)
color_values <-c("#d5231d","#3777ac","#4ea64a","#8e4c99","#e88f18",
                 "#e47faf","#b698c5","#a05528","#58a6d6","#1f2d6f",
                 "#279772","#add387","#d9b71a")

p <- ggtree(Tree, layout = "fan", aes(color = group), 
            branch.length = "none", size = 0.8, open.angle = 10) +
  theme(legend.position = "right") + xlim(2, NA) +
  scale_color_manual(values = color_values,
                     breaks = c(" g__Immundisolibacter"," g__Mesorhizobium"," g__Burkholderia",
                                " g__Methylocystis"," g__Paenibacillus"," g__Pseudacidovorax",
                                " g__Skermanella", " g__Clostridium", " g__Bradyrhizobium",
                                " g__Desulfovibrio", " g__Nostoc"," g__Rhizobium","Others"), 
                     name = "Genus")

###Add circular heat map
p1<-gheatmap(p,logFC, offset=2, 
             width=.3, colnames_angle=95,legend_title = "log2(FC)", font.size = 2) +
  scale_fill_gradient2(low = "#5AC6B4", mid = "white", high = "#F08C8C", 
                       midpoint = 0, name = "log2(FC)")


####################Ecological niche-related calculations##################
##niche width
library(spaa)
library(tidyverse)
library(reshape2)
otu = read.xlsx('Diazotroph.xlsx',1,rowNames = TRUE)
otu <- otu[rowSums(otu != 0) > 0, ]
OTU <-t(otu)

group_labels <- c("CK", "N1", "N5", "W", "WN")
groups <- rep(group_labels, each = 3)  # 定义分组规则


niche_list <- map(seq_along(group_labels), ~ {
  group_data <- OTU[groups == group_labels[.x], ]
  niche_result <- niche.width(group_data, method = "levins") %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("OTUID")
  
  colnames(niche_result)[2] <- group_labels[.x]
  return(niche_result)
})

merged_df <- reduce(niche_list, full_join, by = "OTUID")

df <- merged_df %>% 
  melt(id.vars = "OTUID", value.name = "Value",  variable.name = "Treatment")
df$Value[is.infinite(df$Value)] <- NA

D<-aov(df$Value~df$Treatment,data = df)
anova(D)
P <-LSD.test(D,"df$Treatment")
mar<-P$groups

rownamemar<-row.names(mar)
newmar<-data.frame(rownamemar,mar[1],mar$groups)
sort<-newmar[order(newmar$rownamemar),]
Treatment<-row.names(P$means)
mean<-P$means[,1]
sd<-P$means[,2]
marker<-sort$mar.groups
plotdata<-data.frame(Treatment,mean,sd,marker)

top.mar=0.2
right.mar=0.2
bottom.mar=0.2
left.mar=0.2
mytheme<-theme_bw()+
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        panel.grid=element_blank(),
        panel.border=element_rect(linewidth=1),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))
idcol<-(c("gray","#A7D3D4","#009B9E", "#E4C1D9",  "#C75DAB"))
p<-ggplot(df,aes(x=Treatment,y=Value,color=Treatment))+
  geom_jitter(width=0.15,size=2,show.legend=F)+
  scale_colour_manual(values=idcol)+
  stat_summary(fun = "mean", geom = "point", shape = 8, size = 6) +
  mytheme+
  theme(legend.position = "NA")+
  geom_text(data = plotdata, aes(x = Treatment, y = mean, label = marker), 
            size = 5, color = "black")

####the generalist,opportunist and specialist is identified.
library(EcolUtils)
library(tidyverse)

set.seed(123)

spec_gen <- spec.gen(OTU, niche.width.method = 'levins', perm.method = 'quasiswap', n = 1000, probs = c(0.025, 0.975))
write.csv(spec_gen,'spec_gen.csv')

####################hp-RDA analysis##############
sampledata <- read.xlsx('Diazotroph.xlsx',1,rowNames = TRUE)
Env = read.xlsx('Env.xlsx',1,rowNames = TRUE)
sampledata  <- sampledata [rowSums(sampledata  != 0, na.rm = TRUE) > 0 ,]
spec <- read.xlsx('spec_gen.xlsx',1,rowNames = TRUE)
sampledata$sign <- spec[rownames(sampledata), "sign"]
SPECIALIST <- sampledata[sampledata$sign == "SPECIALIST", -16]
SPECIALIST <- t(SPECIALIST)

#Perform Hellinger transformation on the ASV data
SPECIALIST <- decostand(SPECIALIST,method = "hellinger")

RDA <- rda(SPECIALIST, Env[,c(2:8,10:11)], scale = TRUE)

library(rdacca.hp)
mite.rda.hp <- rdacca.hp(SPECIALIST, Env[,c(2:8,10:11)], method = 'RDA', type = 'R2', scale = FALSE)
mite.rda.hp

set.seed(1)
envfit <- envfit(RDA, Env[,c(2:8,10:11)], permutations  = 999)
r <- as.matrix(envfit$vectors$r)
p <- as.matrix(envfit$vectors$pvals)
env.p <- cbind(r,p)
colnames(env.p) <- c("r2","p-value")
KK <- as.data.frame(env.p)
KK

# 1.Extract the I.perc (%) value and the indicator name from mite.rda.hp
mite_rda_hp_data <- data.frame(
  variable = rownames(mite.rda.hp$Hier.part),
  percentage = mite.rda.hp$Hier.part[, "I.perc(%)"]
)
#2.Extract r2 and p-value from KK
KK_data <- data.frame(
  variable = rownames(KK),
  p_value = KK[, "p-value"]
)

#3. Merge two data frames
data_plot <- merge(mite_rda_hp_data, KK_data, by = "variable")

# 4. Set significance and color conditions
data_plot$significance <- ifelse(data_plot$p_value < 0.05, "Significant", "Not Significant")
data_plot$color <- ifelse(data_plot$p_value < 0.05,"#5AC6B4" , "#F08C8C")

# 5. To maintain the original order of the variables
data_plot$variable <- factor(data_plot$variable, levels = unique(data_plot$variable))

# 6. Plot
QQ<-ggplot(data_plot, aes(x = variable, y = percentage, fill = significance)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(percentage, 2)), vjust = -0.5, size = 3.5) +  # 标注百分比
  coord_polar(start = 1.5) +  # 转换为极坐标
  ylim(-3, max(data_plot$percentage) + 3) +  # 自动设置Y轴范围
  theme_minimal() +
  labs(y = "% Explained Variation", x = "", fill = "Significance") +
  scale_fill_manual(values = c("Significant" = "#5AC6B4", "Not Significant" = "#F08C8C")) +
  theme(
    axis.text.x = element_text(size = 10, hjust = 1),  # x轴标签设置
    panel.background = element_rect(fill = "white")
  )

############community assembly############
tax = read.xlsx('spec_gen.xlsx',1,rowNames = TRUE)
Tax<- merge(otu, tax, by = "row.names")
rownames(Tax) <- Tax[,1]
Tax[,1] <- NULL
specialist<-Tax[Tax$sign == "SPECIALIST",]
opportunitist<-Tax[Tax$sign == "NON SIGNIFICANT",]
generalist<-Tax[Tax$sign == "GENERALIST",]

library(iCAMP)
library(picante)
tree <- read.tree('tree.nwk')
OTU <- t(generalist[,c(1:15)])

abundance=0.05

OTU <- OTU[,colSums(OTU)/sum(OTU)>=(abundance/100)]

spid.check=match.name(cn.list=list(OTU=OTU),tree.list=list(tree=tree))
comm=spid.check$OTU
tree=spid.check$tree

pd <- cophenetic(tree) 

#Community assembly analysis based on the method of Ning et al. (2020)
set.seed(123)
icamp.out <- icamp.big(comm = comm, tree = tree, pd.wd = getwd(), ses.cut = 1.96, rc.cut = 0.95, bin.size.limit = 5, rand = 1000, nworker = 4)

bNRI <- bNRI.cm(comm = comm, dis = pd, weighted = TRUE, rand = 1000, sig.index = 'bNRI')

#####the correlation between environmental variables and ecological processes.

library(vegan)   
library(ecodist) 
library(tidyverse) 

# Step 1: Prepare the data
env_data<-Env[,c(2:8,10,11,17)]
bata_NRI<-bNRI$index
env_scale<-data.frame(scale(env_data,center = F))

sample_pairs <- expand.grid(i = 1:nrow(env_scale), j = 1:nrow(env_scale)) %>%
  filter(i < j) 


# Step 2: Group processes by type

sample_pairs$betaNRI <- apply(sample_pairs, 1, function(idx) {
  bata_NRI[idx["i"], idx["j"]]
})

# Classify the process types based on the βNTI values
sample_pairs <- sample_pairs %>%
  mutate(process_type = case_when(
    abs(betaNRI) > 2 ~ "Deterministic",
    abs(betaNRI) <= 2 ~ "Stochastic"
  ))


# Step 3: Calculate the differences in environment variables
env_vars <- colnames(env_scale)

for(var in env_vars) {
  sample_pairs[[paste0(var, "_diff")]] <- apply(sample_pairs, 1, function(idx) {
    i_index <- as.integer(idx["i"])
    j_index <- as.integer(idx["j"])
    abs(env_scale[i_index, var] - env_scale[j_index, var])
  })
}

# Step 4: Conduct relevant analyses in groups
results <- list()  


for(var in env_vars) {
  # Deterministic process analysis
  deter_data <- sample_pairs %>% 
    filter(process_type == "Deterministic")
  
  if(nrow(deter_data) > 0) {
    cor_deter <- cor.test(
      x = deter_data$betaNRI,
      y = deter_data[[paste0(var, "_diff")]],
      method = "spearman"
    )
  } else {
    cor_deter <- list(estimate = NA, p.value = NA)
  }
  
  # Stochastic Process Analysis
  stoch_data <- sample_pairs %>% 
    filter(process_type == "Stochastic")
  
  if(nrow(stoch_data) > 0) {
    cor_stoch <- cor.test(
      x = stoch_data$betaNRI,
      y = stoch_data[[paste0(var, "_diff")]],
      method = "spearman"
    )
  } else {
    cor_stoch <- list(estimate = NA, p.value = NA)
  }
  

  results[[var]] <- data.frame(
    Variable = var,
    Process = c("Deterministic", "Stochastic"),
    Rho = c(cor_deter$estimate, cor_stoch$estimate),
    P.value = c(cor_deter$p.value, cor_stoch$p.value)
  )
}

results_df <- bind_rows(results)

# Step 5: FDR Correction
results_df <- results_df %>%
  mutate(FDR_adjusted = p.adjust(P.value, method = "BH")) %>%
  mutate(Significant = ifelse(FDR_adjusted < 0.05, "*", ""))


# Step 6: Result Visualization
ggplot(results_df, aes(x = Variable, y = Rho, fill = Process)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = Significant, y = Rho), 
            position = position_dodge(0.9)) +
  scale_fill_manual(values = c( "Deterministic" = "#F08C8C", "Stochastic" = "#5AC6B4")) +
  labs(title = "Correlation between environmental variables and ecological processes",
       x = " environmental variables",
       y = "Spearman's Rho") +
  theme_minimal()+
  theme( panel.grid.major = element_blank() )


##############SEM##############
library(piecewiseSEM)
dat = read.xlsx('Env.xlsx',2,rowNames = TRUE)
all<-data.frame(scale(dat,center = F))

####NL/NH/W
model <- psem(
  lm(ANFR~pH+AP+Biomass,data=all),
  lm(Biomass~NH+W,data=all),
  lm(pH~NH,data=all),
  lm(AP~NL,data=all))
summary(model)

dat = read.xlsx('Env.xlsx',3,rowNames = TRUE)
all<-data.frame(scale(dat,center = F))

#####WN
model <- psem(
  lm(Biomass~N:W,data=all),
  lm(NO3N~Biomass,data=all),
  lm(SM~Biomass,data=all),
  lm(NON.SIGNIFICANT~SPECIALIST+NO3N,data=all),
  lm(SPECIALIST~Biomass,data=all),
  lm(ANFR~SM+NO3N+N:W+NON.SIGNIFICANT+SPECIALIST,data=all))
summary(model)

#!/usr/bin/env R

library(PathoStat)
library(phyloseq)
library(edgeR)
library(gridExtra)
library(grid)
library(ggplot2)
library("data.table"); packageVersion("data.table")

load("/Users/mariangela/Documents/Microbiome/otu_table_exp1.RData") #ps_Susan2018
load("/Users/mariangela/Documents/Microbiome/otu_table_exp2.RData") #psRDP0

ls()
#[1] "ps_Susan2018" "psRDP0"

ps_Susan2018
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10849 taxa and 45 samples ]
#sample_data() Sample Data:       [ 45 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 10849 taxa by 6 taxonomic ranks ]
#ntaxa(ps_Susan2018)
#[1] 10849
#nsamples(ps_Susan2018)
#45
#sample_names(ps_Susan2018)
#rank_names(ps_Susan2018)
#sample_variables(ps_Susan2018)
#[1] "original_name" "variety"       "symtomatology" "time"     
# 01-enrichment1	Unknown			Unknown 		Unknown
#otu_table(ps_Susan2018)[1:5, 1:5]
#tax_table(ps_Susan2018)[1:5, 1:5]
#taxa_names(ps_Susan2018)[1:2]
#[1] "GCAGCAGTGGGGAATAT..."
#sample_variables(ps_Susan2018)
#[1] "original_name" "variety" "symtomatology" "time"  
#get_variable(ps_Susan2018, "symtomatology")
#sample_data(ps_Susan2018)$Sym= factor(get_variable(ps_Susan2018, "symtomatology") %in% c("Feces", "Mock", "Skin", "Tongue") ))



### Filter out Class and Kindom contaminants, _i.e._ chloroplast, mitochondria

ps <- subset_taxa(ps_Susan2018, Class != "Chloroplast")
ps <- subset_taxa(ps, Kingdom != "Mitochondria")
ps <- subset_taxa(ps, Kingdom != "Chloroplast")

ps
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 9458 taxa and 45 samples ]
#sample_data() Sample Data:       [ 45 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 9458 taxa by 6 taxonomic ranks ]

sample_data(ps)$variety <- NULL
sample_data(ps)$symtomatology <- NULL
sample_data(ps)$time  <- NULL

# REPLACE the old group factor with a new one
newtypes = c("enrichment", "enrichment", "enrichment", "kosakonia", "kosakonia", "kosakonia", "saline", "saline", "saline","enrichment", "enrichment", "enrichment", "kosakonia", "kosakonia", "kosakonia", "saline", "saline", "saline", "enrichment", "enrichment", "enrichment", "kosakonia", "kosakonia", "kosakonia", "saline", "saline", "saline","enrichment", "enrichment", "enrichment", "kosakonia", "kosakonia", "kosakonia", "saline", "saline", "saline", "inoculum", "inoculum", "inoculum", "seed", "seed", "seed", "soil", "soil", "soil" )
sample_data(ps)$group = factor(newtypes)
get_variable(ps, "group")

# REPLACE the old generation factor with a new one
newtime=c("passage_1","passage_1","passage_1","passage_1","passage_1","passage_1","passage_1","passage_1","passage_1","passage_2","passage_2","passage_2","passage_2","passage_2","passage_2","passage_2","passage_2","passage_2","passage_3","passage_3","passage_3","passage_3","passage_3","passage_3","passage_3","passage_3","passage_3","passage_4","passage_4","passage_4","passage_4","passage_4","passage_4","passage_4","passage_4","passage_4","passage_1","passage_1","passage_1","passage_1","passage_1","passage_1","passage_1","passage_1","passage_1")

sample_data(ps)$generation = factor(newtime)
get_variable(ps, "generation")
type = c("generation","generation","generation","generation","generation","generation","generation",
         "generation","generation","generation","generation","generation","generation","generation",
         "generation","generation","generation","generation","generation","generation","generation",
         "generation","generation","generation","generation","generation","generation","generation",
         "generation","generation","generation","generation","generation","generation","generation",
         "generation","cc","cc","cc","cc","cc","cc","cc","cc","cc")
sample_data(ps)$type = factor(type)
get_variable(ps, "type")


### Filter out samples for enrichment group
ps <- subset_samples(ps, group!="enrichment")


#remove singletons 
ps.1 <- prune_taxa(taxa_sums(ps) > 1, ps)
ps.1 <- prune_samples(sample_sums(ps) > 0, ps.1)
ps.1
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 8214 taxa and 45 samples ]
#sample_data() Sample Data:       [ 45 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 8214 taxa by 6 taxonomic ranks ]

#Graphical Summary: Phylum_Composition
png("1_Phylum_Composition_psSusan2018.png", height=20, width=30, res=600, units='cm')
title = "Phylum Composition"
plot_bar(tax_glom(ps.1, taxrank="Phylum"), fill = "Phylum") 
dev.off()



#Prevalence
ps.1prev <- transform_sample_counts(ps.1, function(x) ifelse(x>0, 1, 0))

#Relative Abundance
ps.1relabund <- transform_sample_counts(ps.1, function(x) x/sum(x))



#Kosakonia
ps_K = subset_samples(ps.1, group=="kosakonia")
ps_K <- transform_sample_counts(ps_K, function(x) x/sum(x))
ps_K_abu <-  filter_taxa(ps_K, function(x) sum(x) > .01, TRUE)


kos_glom <- tax_glom(ps_K_abu, taxrank="Genus")
test <- merge_samples(kos_glom, "generation")
test <- transform_sample_counts(test, function(x) x/sum(x))

png("RelativeAbundance_Kosakonia_01_run1_may30.png", height=20, width=30, res=600, units='cm') 
plot_bar(test, y="Abundance", fill="Genus", title="Kosakonia: Relative Abundance > .01") +
 geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
dev.off()

kos_df <- psmelt(kos_glom)
write.table(kos_df, file="kosakonia_abundance_run1.csv", sep = "\t", quote = F, row.names = F, col.names = T)

kos_glom <- tax_glom(ps_K_abu, taxrank="Phylum")
kos_df <- psmelt(kos_glom)
write.table(kos_df, file="kosakonia_abundance_Phylum_run1.csv", sep = "\t", quote = F, row.names = F, col.names = T)


ps_K_abu <-  filter_taxa(ps_K, function(x) sum(x) > .05, TRUE)
png("2_RelativeAbundance_Kosakonia_05.png", height=20, width=30, res=600, units='cm')
plot_bar(ps_K_abu, x="generation", y="Abundance", fill="Genus", title="Kosakonia: Relative Abundance > .05")
dev.off()

#soil
ps_soil = subset_samples(ps.1, group=="soil")
ps_soil <- transform_sample_counts(ps_soil, function(x) x/sum(x))

ps_soil_abu <-  filter_taxa(ps_soil, function(x) sum(x) > .01, TRUE)

soil_glom <- tax_glom(ps_soil_abu, taxrank="Genus")
test <- merge_samples(soil_glom, "group")
test <- transform_sample_counts(test, function(x) x/sum(x))
png("RelativeAbundance_soil_01_run1_may30.png", height=20, width=30, res=600, units='cm')
plot_bar(test, y="Abundance", fill="Genus", title="Soil: Relative Abundance > .01") +
 geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
dev.off()


soil_df <- psmelt(soil_glom)
write.table(soil_df, file="soil_abundance_run1.csv", sep = "\t", quote = F, row.names = F, col.names = T)

png("RelativeAbundance_soil_01_run1.png", height=20, width=30, res=600, units='cm')
plot_bar(tax_glom(ps_soil_abu, taxrank="Genus"), x="generation", y="Abundance", fill="Genus", title="Soil: Relative Abundance > .0005")
dev.off()

ps_soil_abu <-  filter_taxa(ps_soil, function(x) sum(x) > .01, TRUE)
png("2_RelativeAbundance_soil_05.png", height=20, width=30, res=600, units='cm')
plot_bar(ps_Srelabund1, x="original_name", y="Abundance", fill="Genus", title="Soil: Relative Abundance > .05")
dev.off()

#saline
ps_sal = subset_samples(ps.1, group=="saline")
ps_sal <- transform_sample_counts(ps_sal, function(x) x/sum(x))
ps_sal_abu <-  filter_taxa(ps_sal, function(x) sum(x) > .01, TRUE)


sal_glom <- tax_glom(ps_sal_abu, taxrank="Genus")
test <- merge_samples(sal_glom, "generation")
test <- transform_sample_counts(test, function(x) x/sum(x))
png("RelativeAbundance_saline_01_run1_may30.png", height=20, width=30, res=600, units='cm')
plot_bar(test, y="Abundance", fill="Genus", title="Saline: Relative Abundance > .01") +
 geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
dev.off()


sal_df <- psmelt(sal_glom)
write.table(sal_df, file="saline_abundance_run1.csv", sep = "\t", quote = F, row.names = F, col.names = T)

sal_glom <- tax_glom(ps_sal_abu, taxrank="Phylum")
sal_df <- psmelt(sal_glom)
write.table(sal_df, file="saline_abundance_Phylum_run1.csv", sep = "\t", quote = F, row.names = F, col.names = T)


ps_sal_abu <-  filter_taxa(ps_sal, function(x) sum(x) > .05, TRUE)
png("2_RelativeAbundance_saline_05.png", height=20, width=30, res=600, units='cm')
plot_bar(ps_sal_abu, x="generation", y="Abundance", fill="Genus", title="Saline: Relative Abundance > .05")
dev.off()


#inoculum
ps_ino = subset_samples(ps.1, group=="inoculum")
ps_ino <- transform_sample_counts(ps_ino, function(x) x/sum(x))
ps_ino_abu <-  filter_taxa(ps_ino, function(x) sum(x) > .01, TRUE)
levels(sample_data(ps_ino_abu)$generation)
sample_data(ps_ino_abu)$generation = factor(sample_data(ps_ino_abu)$generation, levels = c("first","second","third","fourth"))

ino_glom <- tax_glom(ps_ino_abu, taxrank="Genus")
test <- merge_samples(ino_glom, "group")
test <- transform_sample_counts(test, function(x) x/sum(x))
png("RelativeAbundance_inoculum_01_run1.png", height=20, width=30, res=600, units='cm')
plot_bar(test, y="Abundance", fill="Genus", title="Inoculum: Relative Abundance > .01") +
 geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
dev.off()

ino_df <- psmelt(ino_glom)
write.table(ino_df, file="inoculum_abundance_run1.csv", sep = "\t", quote = F, row.names = F, col.names = T)

ino_glom <- tax_glom(ps_ino_abu, taxrank="Phylum")
ino_df <- psmelt(ino_glom)
write.table(ino_df, file="inoculum_abundance_Phylum_run1.csv", sep = "\t", quote = F, row.names = F, col.names = T)

png("RelativeAbundance_inoculum_01_run1.png", height=20, width=30, res=600, units='cm')
plot_bar(tax_glom(ps_ino_abu, taxrank="Genus"), x="generation", y="Abundance", fill="Genus", title="Inoculum: Relative Abundance > .01")
dev.off()


ps_ino_abu <-  filter_taxa(ps_ino, function(x) sum(x) > .05, TRUE)
png("2_RelativeAbundance_inoculum_05.png", height=20, width=30, res=600, units='cm')
plot_bar(ps_ino_abu, x="original_name", y="Abundance", fill="Genus", title="Inoculum: Relative Abundance > .05")
dev.off()


#seed
ps_seed = subset_samples(ps.1, group=="seed")
ps_seed  <- transform_sample_counts(ps_seed, function(x) x/sum(x))
ps_seed_abu <-  filter_taxa(ps_seed, function(x) sum(x) > .01, TRUE)
levels(sample_data(ps_seed_abu)$generation)
sample_data(ps_seed_abu)$generation = factor(sample_data(ps_seed_abu)$generation, levels = c("first","second","third","fourth"))

seed_glom <- tax_glom(ps_seed_abu, taxrank="Genus")
test <- merge_samples(seed_glom, "group")
test <- transform_sample_counts(test, function(x) x/sum(x))
png("RelativeAbundance_seed_01_run1.png", height=20, width=30, res=600, units='cm')
plot_bar(test, y="Abundance", fill="Genus", title="Seeds: Relative Abundance > .01") +
 geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
dev.off()


png("RelativeAbundance_seed_01_run1.png", height=20, width=30, res=600, units='cm')
plot_bar(tax_glom(ps_seed_abu, taxrank="Genus"), x="original_name", y="Abundance", fill="Genus", title="Seed: Relative Abundance > .01")
dev.off()

ps_seed_abu <-  filter_taxa(ps_seed, function(x) sum(x) > .05, TRUE)
png("2_RelativeAbundance_seed_05.png", height=20, width=30, res=600, units='cm')
plot_bar(ps_seed_abu, x="original_name", y="Abundance", fill="Genus", title="Seed: Relative Abundance > .05")
dev.off()


#PcoA, bray’s distances
logt  = transform_sample_counts(ps.1, function(x) log(1 + x) )
out.pcoa.logt <- ordinate(logt, method = "PCoA", distance = "bray")
evals <- out.pcoa.logt$values$Eigenvalues


png("PcoA_run1.png", height=20, width=30, res=600, units='cm')
plot_ordination(logt, out.pcoa.logt, type = "samples", 
                color = "generation", shape = "group") + labs(col = "generation") +
  coord_fixed(sqrt(evals[2] / evals[1]))
dev.off()

png("Ordination_run1.png", height=20, width=30, res=600, units='cm')
plot_ordination(logt, out.pcoa.logt, type = "species", color = "Phylum") 
coord_fixed(sqrt(evals[2] / evals[1]))
dev.off()


#Cluster dendogram
png("ClusterDendogram_run1.png", height=20, width=30, res=600, units='cm')
hell.tip.labels <- as(get_variable(ps.1, "group"), "character")
# This is the actual hierarchical clustering call, specifying average-linkage clustering
d <- distance(ps.1, method="bray", type="samples")
hell.hclust     <- hclust(d, method="average")
plot(hell.hclust)
dev.off()

#remove soil, seed and inoculum
ps.clean <- subset_samples(ps.1, group!="soil")
ps.clean <- subset_samples(ps.clean, group!="seed")
ps.clean <- subset_samples(ps.clean, group!="inoculum")


#get sample sums: each sample's sequencing depth
sample_data(ps.clean)$sequencing_depth <- sample_sums(ps.clean)
MAP <- data.frame(sample_data(ps.clean))
write.table(MAP, file="sequencing_depth_run1_27may.csv", sep = "\t", quote = F, row.names = F, col.names = T)
OTU <- data.frame(as(otu_table(ps.1), 'matrix'))



#Differential Abundances 
#library("edgeR")
m = t(as(otu_table(ps.clean), "matrix"))
m = m + 1
taxonomy = tax_table(ps.clean, errorIfNULL=FALSE)
if( !is.null(taxonomy) ){
  taxonomy = data.frame(as(taxonomy, "matrix"))
} 
d = DGEList(counts=m, genes=taxonomy, remove.zeros = TRUE)

d = calcNormFactors(d, method="RLE")

#check for division by zero inside `calcNormFactors`
if( !all(is.finite(d$samples$norm.factors)) ){
  stop("Something wrong with edgeR::calcNormFactors on this data,
       non-finite $norm.factors, consider changing `method` argument")
}

d$samples

png("MDS_plot_run1_clean.png", height=50, width=50, res=600, units='cm')
colors = c("black", "green", "blue", "red")
plotMDS(z, labels = sample_data(ps.clean)$original_name,  col = as.numeric(factor(sample_data(ps.clean)$generation)),cex = 2, pch = 19)
legend("center", legend = c("generation 1","generation 2", "generation 3", "generation 4" ), col = colors, pch = 19, ncol = 2, cex=2)
dev.off()


#design 
group = as.factor(sample_data(ps.clean)$group)
#generation = as.factor(sample_data(ps.clean)$generation)
generation = as.factor(c("1","2","3","4"))

df <- data.frame(group, generation)
group <- paste0(df$group, ".", df$generation)
cbind(df,group)

df$group <- relevel(df$group, ref="saline")

design <- model.matrix(~group + generation + group:generation, data=df)
#d = estimateDisp(d, design)
d <- estimateGLMCommonDisp(d,design) 
d <- estimateGLMTrendedDisp(d,design) 
d <- estimateGLMTagwiseDisp(d,design)

fit <- glmQLFit(d, design)

colnames(fit)
#[1] "(Intercept)"                "groupkosakonia"            
#[3] "generation2"                "generation3"               
#[5] "generation4"                "groupkosakonia:generation2"
#[7] "groupkosakonia:generation3" "groupkosakonia:generation4"


#the baseline kosakonia vs saline comparison at generation 1
qlf_g1 <- glmQLFTest(fit, coef=2)
dim(qlf_g1)

# plot
tt_g1 <- topTags(qlf_g1, adjust.method="BH", sort.by="PValue")
res_g1 = tt_g1@.Data[[1]]
alpha = 0.05
sigtab_g1 = res_g1[(res_g1$PValue<= alpha), ]
sigtab_g1 = cbind(as(sigtab_g1, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab_g1), ], "matrix"))
dim(sigtab_g1)

png("kosakonia_vs_saline_gen1_run1.png", height=20, width=30, res=600, units='cm')
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab_g1, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
  ggtitle("Log Fold Change of Significant OTUs in Kosakonia vs Saline at generation 1")
dev.off()


#the baseline kosakonia vs saline comparison at generation 2
qlf_g2 <- glmQLFTest(fit, coef=3)
dim(qlf_g2)
topTags(qlf_g2)

png("kosakonia_vs_saline_gen2_run1.png", height=20, width=30, res=600, units='cm')
tt_g2 <- topTags(qlf_g2, adjust.method="BH", sort.by="PValue")
res_g2 = tt_g2@.Data[[1]]
alpha = 0.05
sigtab_g2 = res_g2[(res_g2$PValue<= alpha), ]
sigtab_g2 = cbind(as(sigtab_g2, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab_g2), ], "matrix"))
dim(sigtab_g2)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab_g2, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Kosakonia vs Saline at generation 2")
dev.off()

#the baseline kosakonia vs saline comparison at generation 3
qlf_g3 <- glmQLFTest(fit, coef=4)
dim(qlf_g3)
topTags(qlf_g3)

png("kosakonia_vs_saline_gen3_run1.png", height=20, width=30, res=600, units='cm')
tt_g3 <- topTags(qlf_g3, adjust.method="BH", sort.by="PValue")
res_g3 = tt_g3@.Data[[1]]
alpha = 0.05
sigtab_g3 = res_g3[(res_g3$PValue<= alpha), ]
sigtab_g3 = cbind(as(sigtab_g3, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab_g3), ], "matrix"))
dim(sigtab_g3)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab_g3, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Kosakonia vs Saline at generation 3")
dev.off()



#the baseline kosakonia vs saline comparison at generation 4
qlf_g4 <- glmQLFTest(fit, coef=5)
dim(qlf_g4)
topTags(qlf_g4)

png("kosakonia_vs_saline_gen4_run1.png", height=20, width=30, res=600, units='cm')
tt_g4 <- topTags(qlf_g4, adjust.method="BH", sort.by="PValue")
res_g4 = tt_g4@.Data[[1]]
alpha = 0.05
sigtab_g4 = res_g4[(res_g4$PValue<= alpha), ]
sigtab_g4 = cbind(as(sigtab_g4, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab_g4), ], "matrix"))
dim(sigtab_g4)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab_g4, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Kosakonia vs Saline at generation 4")
dev.off()


#The last two coefficients, differently kosakonia vs saline, at either of the times. 
qlf_any <- glmQLFTest(fit, coef=c(6,7,8))
dim(qlf_any)
topTags(qlf_any)


png("kosakonia_vs_saline_Any_run1.png", height=20, width=30, res=600, units='cm')
tt_any <- topTags(qlf_any, adjust.method="BH", sort.by="PValue")
res_any = tt_any@.Data[[1]]
alpha = 0.05
sigtab_any = res_any[(res_any$PValue<= alpha), ]
sigtab_any = cbind(as(sigtab_any, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab_any), ], "matrix"))
dim(sigtab_any)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab_any, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Kosakonia vs Saline groups")
dev.off()


group = as.factor(sample_data(ps.clean)$group)
#generation = as.factor(sample_data(ps.clean)$generation)
generation = as.factor(c("1","2","3","4"))

df <- data.frame(group, generation)
group <- paste0(df$group, ".", df$generation)
cbind(df,group)


design2 <- model.matrix(~0+group)

d <- estimateGLMCommonDisp(d,design2) 
d <- estimateGLMTrendedDisp(d,design2) 
d <- estimateGLMTagwiseDisp(d,design2)


###quasi-likelihood F-tests:
fit2 <- glmQLFit(d,design2)
colnames(fit2)
[1] "groupkosakonia.1" "groupkosakonia.2" "groupkosakonia.3" "groupkosakonia.4"
[5] "groupsaline.1"    "groupsaline.2"    "groupsaline.3"    "groupsaline.4" 

contr_kos <- makeContrasts(
         Kos.2vs1 = groupkosakonia.2-groupkosakonia.1,
         Kos.3vs1 = groupkosakonia.3-groupkosakonia.1,
         Kos.4vs1 = groupkosakonia.4-groupkosakonia.1,
         Sal.2vs1 = groupsaline.2-groupsaline.1,
         Sal.3vs1 = groupsaline.3-groupsaline.1,
         Sal.4vs1 = groupsaline.4-groupsaline.1,         
         levels=design2)
         
         
qlf <- glmQLFTest(fit2, contrast=contr_kos[,"Kos.2vs1"])
dim(qlf)

png("kosakonia_2_vs_1_run1.png", height=20, width=30, res=600, units='cm')
tt <- topTags(qlf, adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$PValue<= alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab), ], "matrix"))
dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Kosakonia generation 2 vs 1")
dev.off()



qlf <- glmQLFTest(fit2, contrast=contr_kos[,"Kos.3vs1"])
dim(qlf)

png("kosakonia_3_vs_1_run1.png", height=20, width=30, res=600, units='cm')
tt <- topTags(qlf, adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$PValue<= alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab), ], "matrix"))
dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Kosakonia generation 3 vs 1")
dev.off()



qlf <- glmQLFTest(fit, contrast=contr_kos[,"Kos.4vs1"])
dim(qlf)

png("kosakonia_4_vs_1_run1.png", height=20, width=30, res=600, units='cm')

tt <- topTags(qlf, adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$PValue<= alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab), ], "matrix"))
dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Kosakonia generation 4 vs 1")
dev.off()


######### Kosakonia all

### saline all

group = as.factor(sample_data(ps.clean)$group)
#generation = as.factor(sample_data(ps.clean)$generation)
generation = as.factor(c("1","2","3","4"))

df <- data.frame(group, generation)
group <- paste0(df$group, ".", df$generation)
cbind(df,group)


design2 <- model.matrix(~0+group)
#rownames(design) = sample_data(ps.clean)$original_name
#colnames(design) <- levels(group)

d <- estimateGLMCommonDisp(d,design2) 
d <- estimateGLMTrendedDisp(d,design2) 
d <- estimateGLMTagwiseDisp(d,design2)


#To perform quasi-likelihood F-tests:
fit2 <- glmQLFit(d,design2)
colnames(fit2)

contr_kos <- makeContrasts(
         Kos.2vs1 = groupkosakonia.2-groupkosakonia.1,
         Kos.3vs1 = groupkosakonia.3-groupkosakonia.1,
         Kos.4vs1 = groupkosakonia.4-groupkosakonia.1,
         levels=design2)
         
         
qlf <- glmQLFTest(fit2, contrast=contr_kos)
dim(qlf)

png("kosakonia_over_time_run1.png", height=20, width=30, res=600, units='cm')
tt <- topTags(qlf, adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$PValue<= alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab), ], "matrix"))
dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC.Kos.4vs1, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC.Kos.4vs1, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC.Kos.4vs1, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Saline over time - run1")
dev.off()





##### Saline

qlf <- glmQLFTest(fit2, contrast=contr_kos[,"Sal.2vs1"])
dim(qlf)

png("saline_2_vs_1_run1.png", height=20, width=30, res=600, units='cm')
tt <- topTags(qlf, adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$PValue<= alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab), ], "matrix"))
dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Saline generation 2 vs 1")
dev.off()


qlf <- glmQLFTest(fit2, contrast=contr_kos[,"Sal.3vs1"])
dim(qlf)

png("saline_3_vs_1_run1.png", height=20, width=30, res=600, units='cm')
tt <- topTags(qlf, adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$PValue<= alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab), ], "matrix"))
dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Saline generation 3 vs 1")
dev.off()


qlf <- glmQLFTest(fit2, contrast=contr_kos[,"Sal.4vs1"])
dim(qlf)

png("saline_4_vs_1_run1.png", height=20, width=30, res=600, units='cm')
tt <- topTags(qlf, adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$PValue<= alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab), ], "matrix"))
dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Saline generation 4 vs 1")
dev.off()


### saline all

group = as.factor(sample_data(ps.clean)$group)
#generation = as.factor(sample_data(ps.clean)$generation)
generation = as.factor(c("1","2","3","4"))

df <- data.frame(group, generation)
group <- paste0(df$group, ".", df$generation)
cbind(df,group)


design2 <- model.matrix(~0+group)
#rownames(design) = sample_data(ps.clean)$original_name
#colnames(design) <- levels(group)

d <- estimateGLMCommonDisp(d,design2) 
d <- estimateGLMTrendedDisp(d,design2) 
d <- estimateGLMTagwiseDisp(d,design2)


#To perform quasi-likelihood F-tests:
fit2 <- glmQLFit(d,design2)
colnames(fit2)

contr_sal <- makeContrasts(
         Sal.2vs1 = groupsaline.2-groupsaline.1,
         Sal.3vs1 = groupsaline.3-groupsaline.1,
         Sal.4vs1 = groupsaline.4-groupsaline.1,
         levels=design2)
         
         
qlf <- glmQLFTest(fit2, contrast=contr_sal)
dim(qlf)

png("saline_over_time_run1.png", height=20, width=30, res=600, units='cm')
tt <- topTags(qlf, adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$PValue<= alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab), ], "matrix"))
dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC.Sal.4vs1, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC.Sal.4vs1, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC.Sal.4vs1, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Saline over time - run1")
dev.off()



#Alpha Diversity
#png("AlphaDiversity_run1_beforeFiltering.png", height=20, width=30, res=600, units='cm')
pAlpha = plot_richness(ps.1,
                       shape = "generation",
                       color = "group",
                       measures = c("Shannon", "InvSimpson"),
                       title = "Alpha Diveristy")
pAlpha + geom_point(size = 2) + theme(legend.position="bottom")
#dev.off()

alphadt = data.table(pAlpha$data)
write.table(alphadt, file="AlphaDiversity_run1.csv", sep = "\t", quote = F, row.names = F, col.names = T)


### Beta Diversity: Distances --> generation

ps.clean_beta = transform_sample_counts(ps.clean, function(x) x / sum(x))
# Calculate distances
DistBC = distance(ps.clean_beta, method = "bray")
ordBC = ordinate(ps.clean_beta, method = "PCoA", distance = DistBC)

sample_data(ps.clean_beta)$generation = factor(sample_data(ps.clean_beta)$generation, levels = c("first","second","third","fourth"))

png("BetaDiversity_generation_run1.png", height=30, width=30, res=600, units='cm')
plot_ordination(ps.clean_beta, ordBC, color = "group") + 
  geom_point(mapping = aes(size =  generation, 
                           shape = factor(group))) +
  ggtitle("PCoA: Bray-Curtis")
dev.off()


### Beta Diversity: Distances --> cc
ps.cc <- subset_samples(ps.1, group!="kosakonia")
ps.cc <- subset_samples(ps.cc, group!="saline")

ps.cc = transform_sample_counts(ps.cc, function(x) x / sum(x))

# Calculate distances
DistBC = distance(ps.cc, method = "bray")
ordBC = ordinate(ps.cc, method = "PCoA", distance = DistBC)

sample_data(ps.cc)$generation = factor(sample_data(ps.cc)$generation, levels = c("first","second","third","fourth"))

png("BetaDiversity_cc_run1.png", height=20, width=20, res=600, units='cm')
plot_ordination(ps.cc, ordBC, color = "group") + 
  geom_point(mapping = aes(shape = factor(group))) +
  ggtitle("PCoA: Bray-Curtis")
dev.off()






###########################################################################################

load("/Users/mariangela/Documents/Microbiome/otu_table_exp2.RData") #psRDP0

psRDP0
psr <- subset_taxa(psRDP0, Class != "Chloroplast")
psr <- subset_taxa(psr, Kingdom != "Mitochondria")
psr <- subset_taxa(psr, Kingdom != "Chloroplast")

psr
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 3714 taxa and 49 samples ]
#sample_data() Sample Data:       [ 49 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 4545 taxa by 6 taxonomic ranks ]

#original_name, variety, symtomatology, Plant_part


sample_data(psr)$variety <- NULL
sample_data(psr)$symtomatology <- NULL
sample_data(psr)$Plant_part  <- NULL


# REPLACE the old group factor with a new one
newtypes = c("enrich","enrich","enrich","saline","saline","saline","nitrogen","nitrogen","nitrogen","enrich","enrich","enrich","saline","saline","saline","nitrogen","nitrogen","nitrogen","enrich","enrich","enrich","saline","saline","saline","nitrogen","nitrogen","nitrogen","enrich","enrich","enrich","saline","saline","saline","nitrogen","nitrogen","nitrogen","seedlings","seedlings","seedlings","endophytic","endophytic","endophytic","seeds","seeds","seeds","soil","soil","soil","Undetermined_S0")
sample_data(psr)$group = factor(newtypes)
get_variable(psr, "group")

# REPLACE the old original_name
newnames <- c("1enrich1","1enrich2","1enrich3","1saline1","1saline2","1saline3","1nitrogen1","1nitrogen2","1nitrogen3","2enrich1","2enrich2","2enrich3","2saline1","2saline2","2saline3","2nitrogen1","2nitrogen2","2nitrogen3","3enrich1","3enrich2","3enrich3","3saline1","3saline2","3saline3","3nitrogen1","3nitrogen2","3nitrogen3","4enrich1","4enrich2","4enrich3","4saline1","4saline2","4saline3","4nitrogen1","4nitrogen2","4nitrogen3","seedlings1","seedlings2","seedlings3","inoculum1","inoculum2","inoculum3","seeds1","seeds2","seeds3","soil1","soil2","soil3","Undetermined_S0")
sample_data(psr)$original_name = factor(newnames)

# REPLACE the old generation factor with a new one
newtime=c("passage_1","passage_1","passage_1","passage_1","passage_1","passage_1","passage_1","passage_1","passage_1","passage_2","passage_2","passage_2","passage_2","passage_2","passage_2","passage_2","passage_2","passage_2","passage_3","passage_3","passage_3","passage_3","passage_3","passage_3","passage_3","passage_3","passage_3","passage_4","passage_4","passage_4","passage_4","passage_4","passage_4","passage_4","passage_4","passage_4","passage_1","passage_1","passage_1","passage_1","passage_1","passage_1","passage_1","passage_1","passage_1","passage_1","passage_1","passage_1","Undetermined_S0")
sample_data(psr)$generation = factor(newtime)
get_variable(psr, "generation")

type = c("generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","generation","cc","cc","cc","cc","cc","cc","cc","cc","cc","cc","cc","cc","Undetermined_S0")
sample_data(psr)$type = factor(type)
get_variable(psr, "type")


### Filter out samples for enrichment group
psr <- subset_samples(psr, group!="nitrogen")
#Undetermined_S0
psr <- subset_samples(psr, group!="Undetermined_S0")


#remove singletons 
psr.1 <- prune_taxa(taxa_sums(psr) > 1, psr)
psr.1 <- prune_samples(sample_sums(psr.1) > 0, psr.1)
psr.1
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 3114 taxa and 36 samples ]
#sample_data() Sample Data:       [ 36 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 3114 taxa by 6 taxonomic ranks ]

#psr.2  <- tax_glom(psr.1, taxrank="Genus") 

#ntaxa(psr.1); ntaxa(psr.2)

#Graphical Summary: Phylum_Composition
png("Phylum_Composition_run2.png", height=20, width=30, res=600, units='cm')
title = "Phylum Composition"
plot_bar(tax_glom(psr.1, taxrank="Phylum") , fill = "Phylum") 
dev.off()


#Prevalence
psr.1prev <- transform_sample_counts(psr.1, function(x) ifelse(x>0, 1, 0))

#Relative Abundance
psr.1relabund <- transform_sample_counts(psr.1, function(x) x/sum(x))

#enrich
ps_enrich = subset_samples(psr.1, group=="enrich")
ps_enrich <- transform_sample_counts(ps_enrich, function(x) x/sum(x))

ps_enrich_abu <-  filter_taxa(ps_enrich, function(x) sum(x) > .01, TRUE)


enrich_glom <- tax_glom(ps_enrich_abu, taxrank="Genus")
test <- merge_samples(enrich_glom, "generation")
test <- transform_sample_counts(test, function(x) x/sum(x))
png("RelativeAbundance_enrich_01_run2_may30.png", height=20, width=30, res=600, units='cm')
plot_bar(test, y="Abundance", fill="Genus", title="Enrichment: Relative Abundance > .01") +
 geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
dev.off()

enrich_df <- psmelt(enrich_glom)
write.table(enrich_df, file="enrich_abundance_run2.csv", sep = "\t", quote = F, row.names = F, col.names = T)

enrich_glom <- tax_glom(ps_enrich_abu, taxrank="Phylum")
enrich_df <- psmelt(enrich_glom)
write.table(enrich_df, file="enrich_abundance_Phylum_run2.csv", sep = "\t", quote = F, row.names = F, col.names = T)




png("RelativeAbundance_enrich_01_run2.png", height=20, width=30, res=600, units='cm')
plot_bar(tax_glom(ps_enrich_abu, taxrank="Genus"), x="generation", y="Abundance", fill="Genus", title="Enrich: Relative Abundance > .01")
dev.off()

ps_enrich_abu <-  filter_taxa(ps_enrich, function(x) sum(x) > .05, TRUE)
png("2_RelativeAbundance_enrich_05.png", height=20, width=30, res=600, units='cm')
plot_bar(ps_enrich_abu, x="generation", y="Abundance", fill="Genus", title="Enrich: Relative Abundance > .05")
dev.off()

#saline
ps_saline = subset_samples(psr.1, group=="saline")
ps_saline <- transform_sample_counts(ps_saline, function(x) x/sum(x))

ps_saline_abu <-  filter_taxa(ps_saline, function(x) sum(x) > .01, TRUE)
#levels(sample_data(ps_saline_abu)$generation)
#sample_data(ps_saline_abu)$generation = factor(sample_data(ps_saline_abu)$generation, levels = c("first","second","third","fourth"))

saline_glom <- tax_glom(ps_saline_abu, taxrank="Genus")

test <- merge_samples(saline_glom, "generation")
test <- transform_sample_counts(test, function(x) x/sum(x))
png("RelativeAbundance_saline_01_run2.png", height=20, width=30, res=600, units='cm')
plot_bar(test, y="Abundance", fill="Genus", title="Saline: Relative Abundance > .01") +
 geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
dev.off()

saline_df <- psmelt(saline_glom)
write.table(saline_df, file="saline_abundance_run2.csv", sep = "\t", quote = F, row.names = F, col.names = T)

saline_glom <- tax_glom(ps_saline_abu, taxrank="Phylum")
saline_df <- psmelt(saline_glom)
write.table(saline_df, file="saline_abundance_Phylum_run2.csv", sep = "\t", quote = F, row.names = F, col.names = T)

png("RelativeAbundance_saline_01_run2.png", height=20, width=30, res=600, units='cm')
plot_bar(tax_glom(ps_saline_abu, taxrank="Genus"), x="generation", y="Abundance", fill="Genus", title="Saline: Relative Abundance > .01")
dev.off()

ps_saline_abu <-  filter_taxa(ps_saline, function(x) sum(x) > .05, TRUE)
png("2_RelativeAbundance_saline_05.png", height=20, width=30, res=600, units='cm')
plot_bar(ps_saline_abu, x="generation", y="Abundance", fill="Genus", title="Saline: Relative Abundance > .05")
dev.off()

#seedlings
ps_seedlings = subset_samples(psr.1, group=="seedlings")
ps_seedlings <- transform_sample_counts(ps_seedlings, function(x) x/sum(x))

ps_seedlings_abu <-  filter_taxa(ps_seedlings, function(x) sum(x) > .01, TRUE)
levels(sample_data(ps_seedlings_abu)$generation)
sample_data(ps_seedlings_abu)$generation = factor(sample_data(ps_seedlings_abu)$generation, levels = c("first","second","third","fourth"))

seedlings_glom <- tax_glom(ps_seedlings_abu, taxrank="Genus")
test <- merge_samples(seedlings_glom, "group")
test <- transform_sample_counts(test, function(x) x/sum(x))
png("RelativeAbundance_seedlings_glom_01_run2.png", height=20, width=30, res=600, units='cm')
plot_bar(test, y="Abundance", fill="Genus", title="Seedlings: Relative Abundance > .01") +
 geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
dev.off()

png("RelativeAbundance_Seedlings_01_run2.png", height=20, width=30, res=600, units='cm')
plot_bar(tax_glom(ps_seedlings_abu, taxrank="Genus"), x="generation", y="Abundance", fill="Genus", title="Seedlings: Relative Abundance > .01")
dev.off()

ps_seedlings_abu <-  filter_taxa(ps_seedlings, function(x) sum(x) > .05, TRUE)
png("2_RelativeAbundance_Seedlings_05.png", height=20, width=30, res=600, units='cm')
plot_bar(ps_seedlings_abu, x="generation", y="Abundance", fill="Genus", title="Seedlings: Relative Abundance > .05")
dev.off()

#endophytic
ps_endophytic = subset_samples(psr.1, group=="endophytic")
ps_endophytic <- transform_sample_counts(ps_endophytic, function(x) x/sum(x))

ps_endophytic_abu <-  filter_taxa(ps_endophytic, function(x) sum(x) > .01, TRUE)

levels(sample_data(ps_endophytic_abu)$generation)
sample_data(ps_endophytic_abu)$generation = factor(sample_data(ps_endophytic_abu)$generation, levels = c("first","second","third","fourth"))

endophy_glom <- tax_glom(ps_endophytic_abu, taxrank="Genus")
test <- merge_samples(endophy_glom, "group")
test <- transform_sample_counts(test, function(x) x/sum(x))
png("RelativeAbundance_endophytic_01_run2.png", height=20, width=30, res=600, units='cm')
plot_bar(test, y="Abundance", fill="Genus", title="Endophytic: Relative Abundance > .01") +
 geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
dev.off()

endophy_df <- psmelt(endophy_glom)
write.table(endophy_df, file="endophytic_abundance_run2.csv", sep = "\t", quote = F, row.names = F, col.names = T)

endophy_glom <- tax_glom(ps_endophytic_abu, taxrank="Phylum")
endophy_df <- psmelt(endophy_glom)
write.table(endophy_df, file="endophytic_abundance_Phylum_run2.csv", sep = "\t", quote = F, row.names = F, col.names = T)


png("RelativeAbundance_Endophytic_01_run2.png", height=20, width=30, res=600, units='cm')
plot_bar(tax_glom(ps_endophytic_abu,taxrank="Genus"), x="generation", y="Abundance", fill="Genus", title="Endophytic: Relative Abundance > .01")
dev.off()

ps_endophytic_abu <-  filter_taxa(ps_endophytic, function(x) sum(x) > .05, TRUE)
png("RelativeAbundance_Endophytic_05.png", height=20, width=30, res=600, units='cm')
plot_bar(tax_glom(ps_endophytic_abu,taxrank="Genus"), x="generation", y="Abundance", fill="Genus", title="Endophytic: Relative Abundance > .05")
dev.off()

#seeds
ps_seeds = subset_samples(psr.1, group=="seeds")
ps_seeds <- transform_sample_counts(ps_seeds, function(x) x/sum(x))

ps_seeds_abu <-  filter_taxa(ps_seeds, function(x) sum(x) > .01, TRUE)
levels(sample_data(ps_seeds_abu)$generation)
sample_data(ps_seeds_abu)$generation = factor(sample_data(ps_seeds_abu)$generation, levels = c("first","second","third","fourth"))

seeds_glom <- tax_glom(ps_seeds_abu, taxrank="Genus")
test <- merge_samples(seeds_glom, "group")
test <- transform_sample_counts(test, function(x) x/sum(x))
png("RelativeAbundance_seeds_01_run2.png", height=20, width=30, res=600, units='cm')
plot_bar(test, y="Abundance", fill="Genus", title="Seeds: Relative Abundance > .01") +
 geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
dev.off()


png("RelativeAbundance_Seeds_01_run2.png", height=20, width=30, res=600, units='cm')
plot_bar(tax_glom(ps_seeds_abu,taxrank="Genus"), x="generation", y="Abundance", fill="Genus", title="Seeds: Relative Abundance > .01")
dev.off()

ps_seeds_abu <-  filter_taxa(ps_seeds, function(x) sum(x) > .03, TRUE)
png("RelativeAbundance_Seeds_03.png", height=20, width=30, res=600, units='cm')
plot_bar(tax_glom(ps_seeds_abu,taxrank="Genus"), x="generation", y="Abundance", fill="Genus", title="Seeds: Relative Abundance > .03")
dev.off()

#soil
ps_soil = subset_samples(psr.1, group=="soil")
ps_soil <- transform_sample_counts(ps_soil, function(x) x/sum(x))

ps_soil_abu <-  filter_taxa(ps_soil, function(x) sum(x) > .01, TRUE)
levels(sample_data(ps_soil_abu)$generation)
sample_data(ps_soil_abu)$generation = factor(sample_data(ps_soil_abu)$generation, levels = c("first","second","third","fourth"))

soil_glom <- tax_glom(ps_soil_abu, taxrank="Genus")
test <- merge_samples(soil_glom, "group")
test <- transform_sample_counts(test, function(x) x/sum(x))
png("RelativeAbundance_soil_01_run2.png", height=20, width=30, res=600, units='cm')
plot_bar(test, y="Abundance", fill="Genus", title="Soil: Relative Abundance > .01") +
 geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
dev.off()


soil_df <- psmelt(soil_glom)
write.table(soil_df, file="soil_abundance_run2.csv", sep = "\t", quote = F, row.names = F, col.names = T)

soil_glom <- tax_glom(ps_soil_abu, taxrank="Phylum")
soil_df <- psmelt(soil_glom)
write.table(soil_df, file="soil_abundance_Phylum_run2.csv", sep = "\t", quote = F, row.names = F, col.names = T)


png("RelativeAbundance_Soil_01_run2.png", height=20, width=30, res=600, units='cm')
plot_bar(tax_glom(ps_soil_abu, taxrank="Genus"), x="generation", y="Abundance", fill="Genus", title="Soil: Relative Abundance > .01") +
 geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") #to remove the dividing lines that separate OTUs
dev.off()

ps_soil_abu <-  filter_taxa(ps_soil, function(x) sum(x) > .02, TRUE)
png("2_RelativeAbundance_Soil_02.png", height=20, width=30, res=600, units='cm')
plot_bar(ps_soil_abu, x="generation", y="Abundance", fill="Genus", title="Soil: Relative Abundance > .02")
dev.off()


#PCoA
logt  = transform_sample_counts(psr.clean, function(x) log(1 + x) )
out.pcoa.logt <- ordinate(logt, method = "PCoA", distance = "bray")
evals <- out.pcoa.logt$values$Eigenvalues
sample_data(logt)$generation = factor(sample_data(logt)$generation, levels = c("first","second","third","fourth"))

png("PcoA_generation_run2.png", height=20, width=30, res=600, units='cm')
plot_ordination(logt, out.pcoa.logt, type = "samples", 
                color = "generation", shape = "group") + labs(col = "generation") +
  coord_fixed(sqrt(evals[2] / evals[1]))
dev.off()

plot_ordination(logt, out.pcoa.logt, type = "species", color = "Genus") 

png("Ordination_run2.png", height=20, width=30, res=600, units='cm')
plot_ordination(logt, out.pcoa.logt, type = "species", color = "Phylum") 
coord_fixed(sqrt(evals[2] / evals[1]))
dev.off()

#Cluster dendogram
png("ClusterDendogram_run2.png", height=20, width=30, res=600, units='cm')
hell.tip.labels <- as(get_variable(psr.1, "group"), "character")
# This is the actual hierarchical clustering call, specifying average-linkage clustering
d <- distance(psr.1, method="bray", type="samples")
hell.hclust     <- hclust(d, method="average")
plot(hell.hclust)
dev.off()

#remove seedlings, endophytic, seeds, soil
psr.clean <- subset_samples(psr.1, group!="seedlings")
psr.clean <- subset_samples(psr.clean, group!="seeds")
psr.clean <- subset_samples(psr.clean, group!="endophytic")
psr.clean <- subset_samples(psr.clean, group!="soil")

psr.clean
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 3114 taxa and 24 samples ]
#sample_data() Sample Data:       [ 24 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 3114 taxa by 6 taxonomic ranks ]


#get sample sums: each sample's sequencing depth
sample_data(psr.clean)$sequencing_depth <- sample_sums(psr.clean)
MAP <- data.frame(sample_data(psr.clean))
write.table(MAP, file="sequencing_depth_run2_27may.csv", sep = "\t", quote = F, row.names = F, col.names = T)
#OTU <- data.frame(as(otu_table(ps.1), 'matrix'))

##get sample sums: each sample's sequencing depth
sample_data(psr.1)$sequencing_depth <- sample_sums(psr.1)
MAP <- data.frame(sample_data(psr.1))
write.table(MAP, file="sequencing_depth_run2_27may.csv", sep = "\t", quote = F, row.names = F, col.names = T)
#OTU <- data.frame(as(otu_table(ps.1), 'matrix'))



#Differential Abundances 
#library("edgeR")
m = t(as(otu_table(psr.clean), "matrix"))
# Add one to protect against overflow, log(0) issues.
m = m + 1
# Define gene annotations (`genes`) as tax_table
taxonomy = tax_table(psr.clean, errorIfNULL=FALSE)
if( !is.null(taxonomy) ){
  taxonomy = data.frame(as(taxonomy, "matrix"))
} 
# Now turn into a DGEList
d = DGEList(counts=m, genes=taxonomy, remove.zeros = TRUE)

# Calculate the normalization factors
d = calcNormFactors(d, method="RLE")
# Check for division by zero inside `calcNormFactors`
if( !all(is.finite(d$samples$norm.factors)) ){
  stop("Something wrong with edgeR::calcNormFactors on this data,
       non-finite $norm.factors, consider changing `method` argument")
}

png("MDS_plot_run2_clean.png", height=50, width=50, res=600, units='cm')
colors = c("black", "green", "blue", "red")
plotMDS(z, labels = sample_data(psr.clean)$original_name,  col = as.numeric(factor(sample_data(psr.clean)$generation)),cex = 2, pch = 19)
legend("center", legend = c("generation 1","generation 2", "generation 3", "generation 4" ), col = colors, pch = 19, ncol = 2, cex=2)
dev.off()



#design 
group = as.factor(sample_data(psr.clean)$group)
#generation = as.factor(sample_data(psr.clean)$generation)
generation = as.factor(c("1","2","3","4"))

df <- data.frame(group, generation)
group <- paste0(df$group, ".", df$generation)
cbind(df,group)

df$group <- relevel(df$group, ref="saline")

design <- model.matrix(~group + generation + group:generation, data=df)
#d = estimateDisp(d, design)
d <- estimateGLMCommonDisp(d,design) 
d <- estimateGLMTrendedDisp(d,design) 
d <- estimateGLMTagwiseDisp(d,design)

fit <- glmQLFit(d, design)

#This formula is primarily useful as a way to conduct an overall test for interaction.  
#Thecoefficient names are:
colnames(fit)
#[1] "(Intercept)"             "groupenrich"            
#[3] "generation2"             "generation3"            
#[5] "generation4"             "groupenrich:generation2"
#[7] "groupenrich:generation3" "groupenrich:generation4"



#the baseline enrich vs saline comparison at generation 1
qlf_g1 <- glmQLFTest(fit, coef=2)
dim(qlf_g1)

# plot
tt_g1 <- topTags(qlf_g1, adjust.method="BH", sort.by="PValue")
res_g1 = tt_g1@.Data[[1]]
alpha = 0.05
sigtab_g1 = res_g1[(res_g1$PValue<= alpha), ]
sigtab_g1 = cbind(as(sigtab_g1, "data.frame"), as(tax_table(psr.clean)[rownames(sigtab_g1), ], "matrix"))
dim(sigtab_g1)

png("enrich_vs_saline_gen1_run2.png", height=20, width=30, res=600, units='cm')
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab_g1, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
  ggtitle("Log Fold Change of Significant OTUs in Enrichment vs Saline at generation 1")
dev.off()


#the baseline kosakonia vs saline comparison at generation 2
qlf_g2 <- glmQLFTest(fit, coef=3)
dim(qlf_g2)
topTags(qlf_g2)

png("enrich_vs_saline_gen2_run2.png", height=20, width=30, res=600, units='cm')
tt_g2 <- topTags(qlf_g2, adjust.method="BH", sort.by="PValue")
res_g2 = tt_g2@.Data[[1]]
alpha = 0.05
sigtab_g2 = res_g2[(res_g2$PValue<= alpha), ]
sigtab_g2 = cbind(as(sigtab_g2, "data.frame"), as(tax_table(psr.clean)[rownames(sigtab_g2), ], "matrix"))
dim(sigtab_g2)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab_g2, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Enrichment vs Saline at generation 2")
dev.off()

#the baseline kosakonia vs saline comparison at generation 3
qlf_g3 <- glmQLFTest(fit, coef=4)
dim(qlf_g3)
topTags(qlf_g3)

png("enrich_vs_saline_gen3_run2.png", height=20, width=30, res=600, units='cm')
tt_g3 <- topTags(qlf_g3, adjust.method="BH", sort.by="PValue")
res_g3 = tt_g3@.Data[[1]]
alpha = 0.05
sigtab_g3 = res_g3[(res_g3$PValue<= alpha), ]
sigtab_g3 = cbind(as(sigtab_g3, "data.frame"), as(tax_table(psr.clean)[rownames(sigtab_g3), ], "matrix"))
dim(sigtab_g3)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab_g3, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Enrichment vs Saline at generation 3")
dev.off()



#the baseline kosakonia vs saline comparison at generation 4
qlf_g4 <- glmQLFTest(fit, coef=5)
dim(qlf_g4)


png("enrich_vs_saline_gen4_run2.png", height=20, width=30, res=600, units='cm')
tt_g4 <- topTags(qlf_g4, adjust.method="BH", sort.by="PValue")
res_g4 = tt_g4@.Data[[1]]
alpha = 0.05
sigtab_g4 = res_g4[(res_g4$PValue<= alpha), ]
sigtab_g4 = cbind(as(sigtab_g4, "data.frame"), as(tax_table(psr.clean)[rownames(sigtab_g4), ], "matrix"))
dim(sigtab_g4)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab_g4, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Enrichment vs Saline at generation 4")
dev.off()

** ? ***
#The last two coefficients, differently enrich vs saline, at either of the times. 
qlf_any <- glmQLFTest(fit, coef=c(6:8))
dim(qlf_any)



png("enrich_vs_saline_Any_run2.png", height=20, width=30, res=600, units='cm')
tt_any <- topTags(qlf_any, adjust.method="BH", sort.by="PValue")
res_any = tt_any@.Data[[1]]
alpha = 0.05
sigtab_any = res_any[(res_any$PValue<= alpha), ]
sigtab_any = cbind(as(sigtab_any, "data.frame"), as(tax_table(psr.clean)[rownames(sigtab_any), ], "matrix"))
dim(sigtab_any)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab_any, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Enrichment vs Saline groups")
dev.off()



####
group = as.factor(sample_data(psr.clean)$group)
#generation = as.factor(sample_data(psr.clean)$generation)
generation = as.factor(c("1","2","3","4"))

df <- data.frame(group, generation)
group <- paste0(df$group, ".", df$generation)
cbind(df,group)


design2 <- model.matrix(~0+group)
#rownames(design) = sample_data(psr.clean)$original_name
#colnames(design) <- levels(group)

d <- estimateGLMCommonDisp(d,design2) 
d <- estimateGLMTrendedDisp(d,design2) 
d <- estimateGLMTagwiseDisp(d,design2)


#To perform quasi-likelihood F-tests:
fit2 <- glmQLFit(d,design2)
colnames(fit2)
[1] "groupenrich.1" "groupenrich.2" "groupenrich.3" "groupenrich.4"
[5] "groupsaline.1" "groupsaline.2" "groupsaline.3" "groupsaline.4"

contr_kos <- makeContrasts(
         Kos.2vs1 = groupenrich.2-groupenrich.1,
         Kos.3vs1 = groupenrich.3-groupenrich.1,
         Kos.4vs1 = groupenrich.4-groupenrich.1,
         Sal.2vs1 = groupsaline.2-groupsaline.1,
         Sal.3vs1 = groupsaline.3-groupsaline.1,
         Sal.4vs1 = groupsaline.4-groupsaline.1,         
         levels=design2)
         
         
qlf <- glmQLFTest(fit2, contrast=contr_kos[,"Kos.2vs1"])
dim(qlf)

png("enrich_2_vs_1_run2.png", height=20, width=30, res=600, units='cm')
tt <- topTags(qlf, adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$PValue<= alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(psr.clean)[rownames(sigtab), ], "matrix"))
dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Enrichment generation 2 vs 1")
dev.off()



qlf <- glmQLFTest(fit2, contrast=contr_kos[,"Kos.3vs1"])
dim(qlf)

png("enrich_3_vs_1_run2.png", height=20, width=30, res=600, units='cm')
tt <- topTags(qlf, adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$PValue<= alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab), ], "matrix"))
dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Enrichment generation 3 vs 1")
dev.off()



qlf <- glmQLFTest(fit, contrast=contr_kos[,"Kos.4vs1"])
dim(qlf)

png("enrich_4_vs_1_run2.png", height=20, width=30, res=600, units='cm')

tt <- topTags(qlf, adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$PValue<= alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab), ], "matrix"))
dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Enrichment generation 4 vs 1")
dev.off()


### Enrichment all

group = as.factor(sample_data(psr.clean)$group)
#generation = as.factor(sample_data(psr.clean)$generation)
generation = as.factor(c("1","2","3","4"))

df <- data.frame(group, generation)
group <- paste0(df$group, ".", df$generation)
cbind(df,group)


design2 <- model.matrix(~0+group)
#rownames(design) = sample_data(psr.clean)$original_name
#colnames(design) <- levels(group)

d <- estimateGLMCommonDisp(d,design2) 
d <- estimateGLMTrendedDisp(d,design2) 
d <- estimateGLMTagwiseDisp(d,design2)


#To perform quasi-likelihood F-tests:
fit2 <- glmQLFit(d,design2)
colnames(fit2)
[1] "groupenrich.1" "groupenrich.2" "groupenrich.3" "groupenrich.4"
[5] "groupsaline.1" "groupsaline.2" "groupsaline.3" "groupsaline.4"

contr_kos <- makeContrasts(
         Kos.2vs1 = groupenrich.2-groupenrich.1,
         Kos.3vs1 = groupenrich.3-groupenrich.1,
         Kos.4vs1 = groupenrich.4-groupenrich.1,
         levels=design2)
         
         
qlf <- glmQLFTest(fit2, contrast=contr_kos)
dim(qlf)

png("enrich_over_time_run2.png", height=20, width=30, res=600, units='cm')
tt <- topTags(qlf, adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$FDR<= alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(psr.clean)[rownames(sigtab), ], "matrix"))
dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC.Kos.4vs1, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC.Kos.4vs1, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC.Kos.4vs1, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Enrichment over time")
dev.off()


##### Saline

qlf <- glmQLFTest(fit2, contrast=contr_kos[,"Sal.2vs1"])
dim(qlf)

png("saline_2_vs_1_run2.png", height=20, width=30, res=600, units='cm')
tt <- topTags(qlf, adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$PValue<= alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab), ], "matrix"))
dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Saline generation 2 vs 1")
dev.off()


qlf <- glmQLFTest(fit2, contrast=contr_kos[,"Sal.3vs1"])
dim(qlf)

png("saline_3_vs_1_run2.png", height=20, width=30, res=600, units='cm')
tt <- topTags(qlf, adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$PValue<= alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab), ], "matrix"))
dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Saline generation 3 vs 1")
dev.off()


qlf <- glmQLFTest(fit2, contrast=contr_kos[,"Sal.4vs1"])
dim(qlf)

png("saline_4_vs_1_run2.png", height=20, width=30, res=600, units='cm')
tt <- topTags(qlf, adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$PValue<= alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab), ], "matrix"))
dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Saline generation 4 vs 1")
dev.off()

### saline  all

group = as.factor(sample_data(psr.clean)$group)
#generation = as.factor(sample_data(psr.clean)$generation)
generation = as.factor(c("1","2","3","4"))

df <- data.frame(group, generation)
group <- paste0(df$group, ".", df$generation)
cbind(df,group)


design2 <- model.matrix(~0+group)
#rownames(design) = sample_data(psr.clean)$original_name
#colnames(design) <- levels(group)

d <- estimateGLMCommonDisp(d,design2) 
d <- estimateGLMTrendedDisp(d,design2) 
d <- estimateGLMTagwiseDisp(d,design2)


#To perform quasi-likelihood F-tests:
fit2 <- glmQLFit(d,design2)
colnames(fit2)
#[1] "groupenrich.1" "groupenrich.2" "groupenrich.3" "groupenrich.4"
#[5] "groupsaline.1" "groupsaline.2" "groupsaline.3" "groupsaline.4"

contr_sal <- makeContrasts(
         Sal.2vs1 = groupsaline.2-groupsaline.1,
         Sal.3vs1 = groupsaline.3-groupsaline.1,
         Sal.4vs1 = groupsaline.4-groupsaline.1,
         levels=design2)
         
         
qlf <- glmQLFTest(fit2, contrast=contr_sal)
dim(qlf)

png("saline_over_time_run2.png", height=20, width=30, res=600, units='cm')
tt <- topTags(qlf, adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$FDR<= alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(psr.clean)[rownames(sigtab), ], "matrix"))
dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC.Sal.4vs1, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC.Sal.4vs1, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC.Sal.4vs1, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  ggtitle("Log Fold Change of Significant OTUs in Saline over time")
dev.off()


#Alpha Diversity
#png("AlphaDiversity_run2_beforeFiltering.png", height=20, width=30, res=600, units='cm')
pAlpha2 = plot_richness(psr.1,
                       shape = "generation",
                       color = "group",
                       measures = c("Shannon", "InvSimpson"),
                       title = "Alpha Diveristy")
#pAlpha + geom_point(size = 2) + theme(legend.position="bottom")
#dev.off()

alphadt2 = data.table(pAlpha2$data)
write.table(alphadt2, file="AlphaDiversity_run2.csv", sep = "\t", quote = F, row.names = F, col.names = T)



### Beta Diversity: Distances --> generation

psr.clean_beta = transform_sample_counts(psr.clean, function(x) x / sum(x))
# Calculate distances
DistBC = distance(psr.clean_beta, method = "bray")
ordBC = ordinate(psr.clean_beta, method = "PCoA", distance = DistBC)

sample_data(psr.clean_beta)$generation = factor(sample_data(psr.clean_beta)$generation, levels = c("first","second","third","fourth"))

png("BetaDiversity_generation_run2.png", height=30, width=30, res=600, units='cm')
plot_ordination(psr.clean_beta, ordBC, color = "group") + 
  geom_point(mapping = aes(size =  generation, 
                           shape = factor(group))) +
  ggtitle("PCoA: Bray-Curtis")
dev.off()


### Beta Diversity: Distances --> cc
ps.cc <- subset_samples(psr.1, group!="enrich")
ps.cc <- subset_samples(ps.cc, group!="saline")

ps.cc = transform_sample_counts(ps.cc, function(x) x / sum(x))

# Calculate distances
DistBC = distance(ps.cc, method = "bray")
ordBC = ordinate(ps.cc, method = "PCoA", distance = DistBC)

sample_data(ps.cc)$generation = factor(sample_data(ps.cc)$generation, levels = c("first","second","third","fourth"))

png("BetaDiversity_cc_run2.png", height=20, width=20, res=600, units='cm')
plot_ordination(ps.cc, ordBC, color = "group") + 
  geom_point(mapping = aes(shape = factor(group))) +
  ggtitle("PCoA: Bray-Curtis")
dev.off()



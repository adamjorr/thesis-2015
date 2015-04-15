library(ggplot2)
library(reshape)
library(grid)

files = c('hogenom_laged.txt','homologene_laged.txt','pantaxa_corraged.txt')
titles = c('HOGENOM','Homologene','Ensembl Pan-taxa Compara')

svg('dist_plots.svg',width = 21)
pushViewport(viewport(layout = grid.layout(1,3)))
for(i in 1:length(files)){
	df <- read.table(files[i], sep="\t", header=TRUE)
	df.m <- melt(df, measure.vars = c('fam_age','len_age'))
	q <- ggplot(df.m) + geom_density(aes(x=value,colour=variable)) + labs(x = "Age (MY)", y = "Density") + scale_colour_discrete(name = "Age", breaks=c("fam_age","len_age"), labels=c("Family Age", "Length Age")) + ggtitle(titles[i])
	#svg(filename = paste(substr(files[i],1,nchar(files[i])-4),'.svg',sep = ''))
	print(q, vp = viewport(layout.pos.row = 1, layout.pos.col = i))
	#print(q)
}
dev.off()

file <- 'pantaxa_corraged.txt'

df <- read.table(file, sep="\t", header=TRUE)
df.m1 <- melt(df, measure.vars = c('fam_age','corr_age'))
q <- ggplot(df.m1) + geom_density(aes(x=value,colour=variable)) + labs(x = "Age (MY)", y = "Density") + scale_colour_discrete(name = "Age", breaks = c("fam_age", "corr_age"), labels = c("Family Age", "Corrected Family Age")) + ggtitle("Pan-taxa Compara Tree Correction")
svg(filename = paste(substr(file,1,nchar(file)-4),'_fam_correction','.svg',sep = ''))
print(q)

df.m2 <- melt(df, measure.vars = c('len_age','corr_len_age'))
q <- ggplot(df.m2) + geom_density(aes(x=value,colour=variable)) + labs(x = "Age (MY)", y = "Density") + scale_colour_discrete(name = "Age", breaks = c("len_age", "corr_len_age"), labels = c("Length Age", "Corrected Length Age")) + ggtitle("Pan-taxa Compara Tree Correction")
svg(filename = paste(substr(file,1,nchar(file)-4),'_len_correction','.svg',sep = ''))
print(q)

file <- 'calibration_out.txt'
df <- read.table(file, sep="\t", header=FALSE)
q <- ggplot(df) + geom_point(aes(x = V1, y = V2)) + labs( x = 'Percent of Nodes Above', y = 'Node Bootstrap' ) + ggtitle("Calibration")
svg(filename = paste(substr(file,1,nchar(file)-4),'.svg',sep = ''))
print(q)

file <- 'EnsIDsAgeCancerMolGen.csv'
df <- read.csv('EnsIDsAgeCancerMolGen.csv', header=TRUE)
df <- df[c('Gene_Symbol','Cancer_Molecular_Genetics','fam_age')]
df <- subset(df, Cancer_Molecular_Genetics == 'Dom' | Cancer_Molecular_Genetics == 'Rec')
q <- ggplot(df) + geom_density(aes(x = fam_age, color = Cancer_Molecular_Genetics)) + labs(x = "Age (MY)", y = "Density") + scale_colour_discrete(name = "Cancer Gene Type", breaks = c("Dom", "Rec"), labels = c("Dominant", "Recessive")) + ggtitle("Family Age in COSMIC Cancer Genes")
p <- ggplot(df) + geom_boxplot(aes(x = Cancer_Molecular_Genetics, y = fam_age, color = Cancer_Molecular_Genetics)) + labs(x = "Cancer Gene Type", y = "Family Age (MY)") + scale_colour_discrete(name = "Cancer Gene Type", breaks = c("Dom", "Rec"), labels = c("Dominant", "Recessive")) + ggtitle("Family Age in COSMIC Cancer Genes")

svg(filename = paste(substr(file,1,nchar(file)-4),'_distribution','.svg',sep = ''))
print(q)
svg(filename = paste(substr(file,1,nchar(file)-4),'_boxplot','.svg',sep = ''))
print(p)

q()

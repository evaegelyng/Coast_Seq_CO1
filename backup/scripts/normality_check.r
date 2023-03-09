## Checking whether richness, abundance and frequency of phyla are normally distributed
### Load table of summary data
summary<-read.table("rich_abun_freq.tsv",header=T,row.names=1)

## Check normality of mean richness using qqplot
pdf("qqplot_rich.pdf")
qqnorm(summary$log10_mean_rich, pch = 1, frame = FALSE)
qqline(summary$log10_mean_rich, col = "steelblue", lwd = 2)
dev.off()

## Check normality of mean relative abundance using qqplot
pdf("qqplot_abun.pdf")
qqnorm(summary$log10_mean_rel_abund, pch = 1, frame = FALSE)
qqline(summary$log10_mean_rel_abund, col = "steelblue", lwd = 2)
dev.off()

## Check normality of mean relative abundance using qqplot
pdf("qqplot_freq.pdf")
qqnorm(summary$sqrt_sites_occur, pch = 1, frame = FALSE)
qqline(summary$sqrt_sites_occur, col = "steelblue", lwd = 2)
dev.off()
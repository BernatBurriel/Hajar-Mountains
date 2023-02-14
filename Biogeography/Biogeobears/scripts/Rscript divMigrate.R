
### Para usar divMigrate necesitamos un archivo genepop.gen o .arq. Lo creamos con radiator y la funcion genomic_converter
### El strata que necesitamos lleva headers (INDIVIDUALS y ??)
library(radiator)
genomic_converter(data = "qc_randsnp_calotrition_092_75ind_15gen_usnp12.vcf", strata = "strata.txt", output = "genepop", verbose = TRUE)

### Generamos las matrices de migracion con diferentes estadisticos (Gst, D de Jost, etc) y las ploteamos

library(diveRsity)
migRes <- divMigrate(infile = "arnoldi_092_75_15_genepop.gen", outfile = "arnoldi_bootstrap_threshold", filter_threshold = 0.05 , stat = "all", plot_network = TRUE, plot_col = "darkblue", boots = FALSE)
groups <- [(1,2,3,3),(5,6,7,8,9),(10,11,12)]
library(qgraph)
# plot the migration network
pops <- c("A1.1","A1.2","A2","A3","B1","B2","B3","B4","B5","O2.1","O2.2","O3")
df1 <- zip(pops,migRes[[3]])


colnames (migRes[[3]]) <- c("A1.1","A1.2","A2","A3","B1","B2","B3","B4","B5","O2.1","O2.2","O3")
nm <- migRes[[3]]
qgraph(migRes[[3]], 
       edge.labels = FALSE,
       threshold=0.05, 
       layout="spring", theme="Hollywood",
       sampleSize = TRUE, edge.color = "navy",
       esize=4, curveAll=TRUE,
       asize=2.7, curve=3,
       border.width=1, vsize=7,
       labels = colnames(migRes))

?qgraph

rm(list = ls())


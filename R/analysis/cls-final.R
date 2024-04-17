
# final clustering

source(paste(path, "Nextcloud/basilica/scripts/cluster-analysis/utils_plot_clusters.R", sep = ""))
source(paste(path, "Nextcloud/basilica/scripts/cluster-analysis/utils_clusters.R", sep = "")) # compute.scores

p3 <- plot.cluster.freq(x = x)

# plot exposure (all clusters)
p4 <- basilica:::plot_exposures(x, types = "SBS")

# plot exposure for specified cluster
# plot.cluster.exposure(x = x, clusterName = "G0", threshold = 0.05)

G0_exp <- plot.cluster.exposure(x = x, clusterName = "G0", threshold = 0.05)
G1_exp <- plot.cluster.exposure(x = x, clusterName = "G1", threshold = 0.05)
G7_exp <- plot.cluster.exposure(x = x, clusterName = "G7", threshold = 0.05)
G10_exp <- plot.cluster.exposure(x = x, clusterName = "G10", threshold = 0.05)
G11_exp <- plot.cluster.exposure(x = x, clusterName = "G11", threshold = 0.05)
G15_exp <- plot.cluster.exposure(x = x, clusterName = "G15", threshold = 0.05)
G19_exp <- plot.cluster.exposure(x = x, clusterName = "G19", threshold = 0.05)


G0_score <- plot.cluster.score(x = x, clusterName = "G0", threshold = 0.05)
G1_score <- plot.cluster.score(x = x, clusterName = "G1", threshold = 0.05)
G7_score <- plot.cluster.score(x = x, clusterName = "G7", threshold = 0.05)
G10_score <- plot.cluster.score(x = x, clusterName = "G10", threshold = 0.05)
G11_score <- plot.cluster.score(x = x, clusterName = "G11", threshold = 0.05)
G15_score <- plot.cluster.score(x = x, clusterName = "G15", threshold = 0.05)
G19_score <- plot.cluster.score(x = x, clusterName = "G19", threshold = 0.05)


p5 <- plot.score.heatmap(x = x, threshold = 0.05)

#===============================================================================

p3
p4

G0_exp
G1_exp
G7_exp
G10_exp
G11_exp
G15_exp
G19_exp

G0_score
G1_score
G7_score
G10_score
G11_score
G15_score
G19_score

p5



map$basilica_serena

basilica:::plot_exposures(x, sort_by = "SBS13")


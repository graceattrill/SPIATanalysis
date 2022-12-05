plot_distance_heatmap_GAscaled<-function(phenotype_distances_result, metric = "mean") 
{
  Reference <- Nearest <- Mean <- Std.Dev <- Median <- NULL
  if (metric == "mean") {
    limit <- c(0,2500)
    g <- ggplot(phenotype_distances_result, aes(x = Reference, 
                                                y = Nearest, fill = Mean)) + geom_tile() + xlab("Cell of interest (COI)") + 
      ylab("Nearest cell to COI") + theme(panel.grid.major = element_blank(), 
                                          panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), 
                                          axis.text.x = element_text(angle = 90, hjust = 1, 
                                                                     vjust = 0.5)) + scale_fill_viridis_c(limits = limit, 
                                                                                                          direction = -1)
    print(g)
  }
  else if (metric == "std.dev") {
    limit <- c(0,3000)
    g <- ggplot(phenotype_distances_result, aes(x = Reference, 
                                                y = Nearest, fill = Std.Dev)) + geom_tile() + xlab("Cell of interest (COI)") + 
      ylab("Nearest cell to COI") + theme(panel.grid.major = element_blank(), 
                                          panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), 
                                          axis.text.x = element_text(angle = 90, hjust = 1, 
                                                                     vjust = 0.5)) + scale_fill_viridis_c(limits = limit, 
                                                                                                          direction = -1)
    print(g)
  }
  else if (metric == "median") {
    limit <- c(4000)
    g <- ggplot(phenotype_distances_result, aes(x = Reference, 
                                                y = Nearest, fill = Median)) + geom_tile() + xlab("Cell of interest (COI)") + 
      ylab("Nearest cell to COI") + theme(panel.grid.major = element_blank(), 
                                          panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), 
                                          axis.text.x = element_text(angle = 90, hjust = 1, 
                                                                     vjust = 0.5)) + scale_fill_viridis_c(limits = limit, 
                                                                                                          direction = -1)
    print(g)
  }
}

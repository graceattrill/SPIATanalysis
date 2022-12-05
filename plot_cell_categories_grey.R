plot_cell_categories<-function (spe_object, categories_of_interest = NULL, colour_vector = NULL, 
          feature_colname = "Cell.Type", cex = 1, layered = FALSE) 
{
  if (feature_colname == "Structure" & is.null(categories_of_interest)) {
    categories_of_interest <- c("Border", "Inside", 
                                "Infiltrated.immune", "Outside", "Stromal.immune", 
                                "Internal.margin", "Internal.margin.immune", 
                                "External.margin", "External.margin.immune")
    colour_vector <- c("black", "pink", "purple", 
                       "yellow", "orange", "lightgreen", 
                       "darkgreen", "lightblue", "blue")
  }
  Cell.X.Position <- Cell.Y.Position <- Category <- NULL
  if (class(sce_object) == "SingleCellExperiment" || 
      class(sce_object) == "SummarizedExperiment") {
    formatted_data <- data.frame(SummarizedExperiment::colData(sce_object))
  }
  else formatted_data <- sce_object
  if (length(categories_of_interest) != length(colour_vector)) {
    stop("The colour vector is not the same length as the phenotypes of interest")
  }
  for (category in categories_of_interest) {
    if (!(category %in% unique(formatted_data[[feature_colname]]))) {
      cat_idx <- match(category, categories_of_interest)
      categories_of_interest <- categories_of_interest[-cat_idx]
      colour_vector <- colour_vector[-cat_idx]
      print(paste(category, "cells were not found and not plotted"), 
            sep = "")
    }
  }
  if (any(!formatted_data[[feature_colname]] %in% categories_of_interest)) {
    formatted_data[!formatted_data[[feature_colname]] %in% 
                     categories_of_interest, ][[feature_colname]] <- "OTHER"
  }
  formatted_data$color <- ""
  for (category in categories_of_interest) {
    idx <- which(categories_of_interest == category)
    formatted_data[formatted_data[[feature_colname]] == category, 
                   ]$color <- colour_vector[idx]
  }
  if (any(formatted_data[[feature_colname]] == "OTHER")) {
    formatted_data[formatted_data[[feature_colname]] == "OTHER", 
                   ]$color <- "lightgrey"
    all_categories <- c(categories_of_interest, "OTHER")
    all_colours <- c(colour_vector, "lightgrey")
  }
  else {
    all_categories <- categories_of_interest
    all_colours <- colour_vector
  }
  if (layered) {
    p <- ggplot(formatted_data, aes_string(x = "Cell.X.Position", 
                                           y = "Cell.Y.Position"))
    for (cat in categories_of_interest) {
      p <- p + geom_point(data = formatted_data[formatted_data[[feature_colname]] == 
                                                  cat, ], aes_string(colour = feature_colname), 
                          size = cex)
    }
    p <- p + guides(alpha = "none") + ggtitle(paste("Plot", 
                                                    attr(sce_object, "name"), feature_colname, 
                                                    sep = " ")) + scale_color_manual(breaks = all_categories, 
                                                                                     values = all_colours)
  }
  else {
    p <- ggplot(formatted_data, aes_string(x = "Cell.X.Position", 
                                           y = "Cell.Y.Position")) +
    geom_point(aes_string(colour = feature_colname), size = cex)
    p <- p + guides(alpha = "none") + ggtitle(paste("Plot", 
                                                    attr(sce_object, "name"), feature_colname, 
                                                    sep = " ")) + scale_color_manual(breaks = all_categories, 
                                                                                     values = all_colours)
  }
  print(p)
}
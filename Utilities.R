heatmap_areas <- function(map_w_values,
                          value,
                          scale_col = NULL,
                          scale = NULL,
                          hardcoded_bins = NULL,
                          title = NULL){
  
  map_w_values$to_plot = value
  
  if(is.null(scale_col)){
    scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 30 
  }
  if(is.null(scale)){
    scale = scale_col[c(3, 8, 12, 15, 19, 23, 26, 30)]
  }
  if(is.null(hardcoded_bins)){
    ggplot(data = map_w_values) +  
      geom_sf(aes(fill = to_plot), 
              alpha = 1,
              color="black") + ggtitle(title) + 
      theme(plot.title = element_text(size = 15),
            axis.title.x = element_blank(), #Remove axis and background grid
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank(),
            plot.margin =  unit(c(0, 0, 0, 0), "inches"),
            legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
            panel.spacing = unit(1, 'lines')) +
      guides(fill=guide_legend(title=NULL, reverse = TRUE, label.position = "right")) + #Remove colorbar title
      binned_scale( #Scaling the color
        aesthetics = "fill",
        scale_name = "gradientn",
        palette = function(x) c(scale),
        labels = function(x){x},
        guide = "colorscale")
  } else {
    ggplot(data = map_w_values) +  
      geom_sf(aes(fill = to_plot), 
              alpha = 1,
              color="black") + ggtitle(title) + 
      theme(plot.title = element_text(size = 15),
            axis.title.x = element_blank(), #Remove axis and background grid
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank(),
            plot.margin =  unit(c(0, 0, 0, 0), "inches"),
            legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
            panel.spacing = unit(1, 'lines')) +
      guides(fill=guide_legend(title=NULL, reverse = TRUE, label.position = "right")) + #Remove colorbar title
      binned_scale( #Scaling the color
        aesthetics = "fill",
        scale_name = "gradientn",
        palette = function(x) c(scale),
        labels = function(x){x},
        breaks = hardcoded_bins,
        guide = "colorscale")
    }
}
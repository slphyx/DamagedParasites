library(ggplot2)


plot.counts <- function(extractedfilename, id = NULL, split = FALSE) {
  tmp <- read.csv(extractedfilename)
  dlim <- min(tmp$obdat)

  if (split == TRUE) {
    outgraph <- 
      ggplot(data = tmp, aes(x = obtime, y = obdat)) +
      geom_ribbon(data = tmp, aes(x = obtime, ymax = above, ymin = below), fill = "pink", alpha = .5) +
      xlab("time (hrs)") +
      ylab("log10 circulating parasites") +
      geom_point(data = tmp, color = "black") +
      geom_line(data = tmp, aes(y = obdat, color = 'Observed')) +
      geom_point(data = tmp, aes(y = median), color = "red") +
      geom_line(data = tmp, aes(y = median, color = "Predicted Median")) +
      geom_point(data = tmp, aes(x = obtime, y = mediansplit), color = "blue") +
      geom_line(data = tmp, aes(x = obtime, y = mediansplit, color = "Splitting dose")) +
      geom_ribbon(data = tmp, aes(x = obtime, ymax = abovesplit, ymin = belowsplit), fill = "lightblue", alpha = .5) +
      scale_color_manual("", breaks = c("Observed", "Predicted Median", "Splitting dose"), 
                         values = c('Observed' = 'black', 'Predicted Median' = "red", "Splitting dose"="blue")) +
      geom_hline(yintercept = dlim, linetype = "longdash", color = "darkgrey") +
      geom_text(aes(10, dlim, label = "detection limit", vjust = 1), size = 3, color = "grey", fontface = "italic") +
      ggtitle(id) +
      theme_bw() +
      theme(legend.position = "bottom", legend.key = element_blank(), plot.title = element_text(hjust = 0.5)
            , plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))+
      coord_cartesian(expand = FALSE,ylim = c(-1,14), xlim = c(-5,115)) 
      coord_fixed(ratio = 3.5)
    
  } else {
    outgraph <- 
      ggplot(data = tmp, aes(x = obtime, y = obdat)) +
      geom_ribbon(data = tmp, aes(x = obtime, ymax = above, ymin = below), fill = "pink", alpha = .5) +
      xlab("time (hrs)") +
      ylab("log10 circulating parasites") +
      geom_point(data = tmp, color = "black") +
      geom_line(data = tmp, aes(y = obdat, color = 'Observed')) +
      geom_point(data = tmp, aes(y = median), color = "red") +
      geom_line(data = tmp, aes(y = median, color = "Predicted Median")) +
      geom_hline(yintercept = dlim, linetype = "longdash", color = "darkgrey") +
      geom_text(aes(10, dlim, label = "detection limit", vjust = 1), size = 3, color = "grey", fontface = "italic") +
      scale_color_manual("", breaks = c("Observed", "Predicted Median"), values = c('Observed' = 'black', 'Predicted Median' = "red")) +
      ggtitle(id) +
      theme_bw() +
      theme(legend.position = "bottom", legend.key = element_blank(), plot.title = element_text(hjust = 0.5)
            , plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))+
      xlim(-5,115)+ylim(-1,14)
  }
  
  outgraph 
}






cumulativeGraph <- function(hci) {
  library(ggplot2)
  d <- data.frame(hci$height, 1:length(hci$height))
  colnames(d) <- c("height", "cumul")
  
  ggplot(d, aes(x = 1-d$height, y = d$cumul)) + geom_line()
  
  plot(d)
}
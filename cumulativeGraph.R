cumulativeGraph <- function(cgraph.d, filen, grid=false) {
  # pass in the dataframe from the global environment
  #cgraph.d <<- data.frame(hci$height, 1:length(hci$height))
  print(cgraph.d)
  #colnames(cgraph.d) <- c("threshold", "cumul")
  print('ok')
  if(!grid){
    p <- ggplot(cgraph.d, aes(x = cgraph.d$height, y = cgraph.d$cumul), environment=environment()) + geom_line()
  } else {
    p <- ggplot(cgraph.d, aes(x = cgraph.d$height, y = cgraph.d$cumul), environment=environment()) + geom_line() + facet_wrap(~region)
  }
  ggsave(filename=filen, plot=p)
  print('ok')
  return(p)
}
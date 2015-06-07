cumulativeGraph <- function(cgraph.d, filen) {
  #print(cgraph.d) # for debugging only
  p <- ggplot(cgraph.d, aes(x = cgraph.d$height, y = cgraph.d$cumul), environment=environment()) + geom_line()
  ggsave(filename=filen, plot=p)
  return(p)
}
showHistPlot <- function(input){
  dist <- rnorm(input$obs)
  return(hist(dist))
}
  
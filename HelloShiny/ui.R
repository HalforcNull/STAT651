library(shiny)
source('customPackages/tryUINest.R')


# Define UI for application that plots random distributions 
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Hello Shiny!"),
  
  # Sidebar with a slider input for number of observations
  sidebarLayout(
    myfun,
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
))

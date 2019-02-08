library(shiny)



ui <- navbarPage("Navbar page", id = "tabs",
				 tabPanel("Home",
				 		 actionButton("hideTab", "Hide 'Foo' tab"),
				 		 actionButton("showTab", "Show 'Foo' tab"),
				 		 actionButton("hideMenu", "Hide 'More' navbarMenu"),
				 		 actionButton("showMenu", "Show 'More' navbarMenu"),
				 		 textOutput('testingTextOutput1'),
				 		 textOutput('testingTextOutput2'),
				 		 textOutput('testingTextOutput3'),
				 		 numericInput("myInput", label = h5("My input"), value = 10)
				 		 
				 ),
				 tabPanel("Foo", "This is the foo tab"),
				 tabPanel("Bar", "This is the bar tab")
)


server <- function(input, output, session) {
	observeEvent(input$hideTab, {
		hideTab(inputId = "tabs", target = "Foo")
	})
	
	observeEvent(input$showTab, {
		showTab(inputId = "tabs", target = "Foo")
	})
	
	observeEvent(input$hideMenu, {
		hideTab(inputId = "tabs", target = "More")
	})
	
	observeEvent(input$showMenu, {
		showTab(inputId = "tabs", target = "More")
	})
	
	output$testingTextOutput1 <-  reactive({renderText(input$myInput)})
	#output$testingTextOutput2 <- input$myInput
	output$testingTextOutput3 <- isolate({renderText(input$myInput)})
}

shinyApp(ui, server)
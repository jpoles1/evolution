source("base.R")
shinyServer(function(input, output) {
  
  output$charts <- renderPlot({ 
    evolve(pop=input$pop, gen=input$gen, p=input$p, wpp=input$wpp, wpq=input$wpq, wqq=input$wqq, pmut=input$pmut, qmut=input$qmut)
  })
  
})
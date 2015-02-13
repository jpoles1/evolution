stepsize = .000001
shinyUI(fluidPage(
  titlePanel("Single Allele Evolution Calculator"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Calculate Genotype and Allele frequency for an allele over time under selection and mutation pressures."),
      sliderInput("pop", label = "Size of Population?", value=100, min=1, max=10000),
      sliderInput("gen", label = "How many generations?", value=10, min=1, max=1000),
      sliderInput("p", label = "Initial frequency of allele P", value=.5, min=0, max=1, step=stepsize),
      sliderInput("wpp", label = "Selection rate for pp", value=1, min=0, max=1, step=stepsize),
      sliderInput("wpq", label = "Selection rate for pq", value=1, min=0, max=1, step=stepsize),
      sliderInput("wqq", label = "Selection rate for qq", value=1, min=0, max=1, step=stepsize),
      sliderInput("pmut", label = "Mutation rate for allele p to q", value=0, min=0, max=1, step=stepsize),
      sliderInput("qmut", label = "Mutation rate for allele q to p", value=0, min=0, max=1, step=stepsize)
    ),
    mainPanel(
      plotOutput("charts")
    )
  )
))
stepsize = .000001
shinyUI(fluidPage(
  titlePanel("Single Allele Evolution Calculator"),
  plotOutput("charts"),
  fluidRow(
    h3("Calculate Genotype and Allele frequency for an allele over time under selection and mutation pressures."),
    column(3, 
      h4("Population Parameters"),
      sliderInput("gen", label = "How many generations?", value=10, min=1, max=1000),
      sliderInput("p", label = "Initial frequency of allele P", value=.5, min=0, max=1, step=stepsize)
    ),
    column(4,
      h4("Selection Rates"),
      sliderInput("wpp", label = "Selection rate for pp", value=1, min=0, max=1, step=stepsize),
      sliderInput("wpq", label = "Selection rate for pq", value=1, min=0, max=1, step=stepsize),
      sliderInput("wqq", label = "Selection rate for qq", value=1, min=0, max=1, step=stepsize)
    ),
    column(4,
      h4("Mutation Rates"),
      sliderInput("pmut", label = "Mutation rate for allele p to q", value=0, min=0, max=1, step=stepsize),
      sliderInput("qmut", label = "Mutation rate for allele q to p", value=0, min=0, max=1, step=stepsize) 
    )
  )
))
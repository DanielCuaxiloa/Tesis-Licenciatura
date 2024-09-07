
library(shiny)
library(shinyjs)
library(shinydashboard)
library(ggplot2)
library(dplyr)

# Cargar los datos de los modelos
Modelo1 <- read.csv(file = "Modelo1.csv", stringsAsFactors = TRUE)
Modelo2 <- read.csv(file = "Modelo2.csv", stringsAsFactors = TRUE)
Modelo3 <- read.csv(file = "Modelo3.csv", stringsAsFactors = TRUE)
Modelo4 <- read.csv(file = "Modelo4.csv", stringsAsFactors = TRUE)
Modelo5 <- read.csv(file = "Modelo5.csv", stringsAsFactors = TRUE)

# Combinar los datos
Datos <- bind_rows(Modelo1, Modelo2, Modelo3, Modelo4, Modelo5) %>% 
  select(-ID1) %>% 
  select(Modelo, Nombre, TestClase1, TestClase2, TestClase3, TestGlobal) %>% 
  rename(
    'TCC General' = TestGlobal,
    'TCC - GTEX Brain' = TestClase1,
    'TCC - TCGA LGG' = TestClase2,
    'TCC - TCGA GM' = TestClase3,
    Model = Nombre
  )

ui <- dashboardPage(
  dashboardHeader(
    title = tags$div(
      tags$img(src = "logoFC85.png", height = "40px", style = "margin-right: 10px;"),
      "Dashboard | Tesis"
    ),
    disable = FALSE
  ),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Menú Principal", tabName = "menu", icon = icon("home")),
      menuItem("Aspectos Teóricos", tabName = "teoricos", icon = icon("book")),
      menuItem("Caso Práctico", tabName = "practico", icon = icon("chart-bar"))
    )
  ),
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .main-content {
          background-color: #f0f0f0;
        }
        .box {
          background-color: #f0f0f0;
          border: 1px solid #d2d6de;
          border-radius: 3px;
          padding: 10px;
        }
      "))
    ),
    useShinyjs(),
    tabItems(
      tabItem(tabName = "menu",
              div(
                class = "centered",
                h2("Menú Principal"),
                p("Selecciona una de las opciones del menú para continuar.")
              )
      ),
      tabItem(tabName = "teoricos",
              tabsetPanel(
                tabPanel("Aprendizaje Supervisado",
                         div(
                           class = "main-content",
                           h3("Aprendizaje Supervisado"),
                           p("...")
                           #...
                         )
                ),
                tabPanel("Algoritmos de Clasificación",
                         div(
                           class = "main-content",
                           h3("Algoritmos de Clasificación"),
                           p("...")
                           #...
                         )
                )
              )
      ),
      tabItem(tabName = "practico",
              tabsetPanel(
                tabPanel("Contexto",
                         div(
                           class = "main-content",
                           h3("Contexto"),
                           p("...")
                           #...
                         )
                ),
                tabPanel("Resultados",
                         div(
                           class = "main-content box",
                           h3("Resultados Modelo"),
                           fluidRow(
                             column(6,
                                    selectInput("ycol", "Métrica", choices = c("TCC General", 
                                                                               "TCC - GTEX Brain", 
                                                                               "TCC - TCGA LGG", 
                                                                               "TCC - TCGA GM"), 
                                                selected = "TCC General")
                             ),
                             column(6,
                                    selectInput("modelos", "Modelo", 
                                                choices = unique(Datos$Modelo), 
                                                multiple = TRUE,
                                                selectize = TRUE) 
                             )
                           ),
                           plotOutput("modelPlot")
                         )
                )
              )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Visualización del modelo
  output$modelPlot <- renderPlot({
    filteredData <- Datos %>% filter(Modelo %in% input$modelos)
    p <- ggplot(filteredData, mapping = aes(x = Model, y = !!sym(input$ycol))) +
      facet_grid(cols = vars(Modelo), scales = "free_x") +
      geom_boxplot(fill = "steelblue3") + 
      theme_bw() + 
      theme(
        panel.background = element_rect(fill = "#f0f0f0", color = NA),
        plot.background = element_rect(fill = "#f0f0f0", color = NA),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"),
        strip.background = element_rect(fill = "#d2d6de"),
        strip.text = element_text(color = "black")
      ) +
      labs(title = paste("Estimación de poder predictivo basado en la métrica", input$ycol, "por Modelo"))
    
    p
  })
  
}

shinyApp(ui = ui, server = server)




library(shiny)
library(shinyjs)
library(shinydashboard)

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
                           p("Aquí se describen los aspectos teóricos del aprendizaje supervisado...")
                           # Añade más contenido teórico sobre aprendizaje supervisado aquí
                         )
                ),
                tabPanel("Algoritmos de Clasificación",
                         div(
                           class = "main-content",
                           h3("Algoritmos de Clasificación"),
                           p("Aquí se describen los aspectos teóricos de los algoritmos de clasificación...")
                           # Añade más contenido teórico sobre algoritmos de clasificación aquí
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
                           p("Aquí se describe el contexto del caso práctico...")
                           # Añade más contenido sobre el contexto del caso práctico aquí
                         )
                ),
                tabPanel("Resultados",
                         div(
                           class = "main-content",
                           h3("Resultados"),
                           p("Aquí se muestran los resultados del caso práctico...")
                           # Añade más contenido sobre los resultados del caso práctico aquí
                         )
                )
              )
      )
    )
  )
)

server <- function(input, output, session) {
  
}

shinyApp(ui = ui, server = server)


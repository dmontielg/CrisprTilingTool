
library(DT)
library(shiny)
library(dplyr)
library(shinycssloaders)
source("get-phase-trancript.R")
source("get-sgRNA.R")
source("locate-base-changes.R")
source("run_locate_modifications.R")


#fixed_cds_path <- "/home/shiny_app/CDS-human.csv"
fixed_cds_path <- "CDS-human.csv"


ui <- fluidPage(
  navbarPage(
    title = ("Crispr Tilling Tool"),

  ),
  sidebarLayout(
    sidebarPanel(
      selectInput("base_change", label = h3("Base change"),choices = c("C->T", "A->G"), multiple = FALSE),
      selectInput("genome", label = ("Genome"),choices = c("hg38"), multiple = FALSE),
      textInput("trans", label = ("Transcript"), value = "ENST00000269305.9"),
      numericInput("flank", label = ("Flank"), value = "30"),
      numericInput("sgrna_len", label = ("sgRNA length"), value = "20", min = 0),
      numericInput("window_start", label = ("Window start index"), value = "4",min = 0),
      numericInput("window_end", label = ("Window end index"), value = "8", min = 1),
      textInput("seq_before", label = ("DNA sequence before sgRNA"), value = "GAGCCTCGTCTCCCACCG"),
      textInput("seq_after", label = ("DNA sequence after sgRNA"), value = "GTTTTGAGACGCATGCTGCA"),
      fluidRow(
        actionButton("submit", label = "Run"),
        downloadButton("downloadData", "download"),
      ),
    ),
    #fluidPage(  fluidRow(DT::dataTableOutput("table"))),
    fluidPage(  fluidRow(shiny::tableOutput("table"))),

  )
)

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 60 * 1024^2)

  CDS_Human <- reactive({
    readr::read_csv(fixed_cds_path)
  })

  observeEvent(input$submit, {
    req(input$trans)
    #showPageSpinner()
    #Sys.sleep(20)
    #hidePageSpinner()
    Genome      <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
    df.target.trans <- CDS_Human() %>%
      dplyr::filter(tx == input$trans)
      if (nrow(df.target.trans) > 0) {
            gene_search <- unique(df.target.trans$gene)
            df_table_out <- run_locate_modifications(df.target.trans,
                                                    gene_search,
                                                    input$trans,
                                                    input$flank,
                                                    input$sgrna_len,
                                                    input$window_start,
                                                    input$window_end,
                                                    input$base_change,
                                                    input$seq_before,
                                                    input$seq_after,
                                                    Genome)
            #print(head(df_table_out))

            #output$table = shiny::renderTable(df_table_out)

            output$table = DT::renderDataTable({df_table_out})

            data_table_react <- reactive({ df_table_out })
            output$downloadData <- downloadHandler(
                    filename = paste("Table", gene_search, "mod", input$base_change,
                                    "flank", input$flank, ".txt", sep = "_", collapse = ""),
                    content = function(file) {
                      write.csv(data_table_react(), file, row.names = F, quote = F)
                    }
                  )
      }
    })
}

shinyApp(ui = ui, server = server)

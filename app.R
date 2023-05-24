

#source functions from function file
# source("functions.R")
source("workflow.R")

# Define UI
ui <- fluidPage(
  shinyFeedback::useShinyFeedback(),
  theme = shinytheme("slate"),
  navbarPage(
    "LUCA Analysis",
    tabPanel(
      "LUCA",
      sidebarPanel(
        tags$h3("Input:"),
        fileInput("edds", "Upload EDDS .csv", accept = ".csv"),
        fileInput(
          "flowjo",
          "Upload flowjo files",
          accept = c('.xls', ".xlsx"),
          multiple = TRUE
        ),
        
        fileInput(
          "dosing",
          "Upload dosing solution readout",
          accept = c('.xls', ".xlsx"),
          multiple = TRUE
        ),
        #textInput('mfi_var','Enter mfi variable',value = 'Geometric Mean : pHAb-A'),
        # verbatimTextOutput('mfi_choices')
        #selectInput('mfi_choices','Select desired MFI variable',choices = ''),
        # actionButton('submit','Submit'),
        textInput(inputId = 'mfi_choices', 'Select variable for MFI', value = "Geometric Mean : pHAb-A"),
        verbatimTextOutput('mfi_choices_list'),
        textInput(
          'control_mabs',
          'Select control antibodies in ascending order, separated by comma'
        ),
        verbatimTextOutput('control_choices'),
        verbatimTextOutput('misc'),
        downloadButton("download_csv", "Download EDDS.csv"), downloadButton("download_pzfx", "Download Graphpad friendly .csv")
      ),
      # sidebarPanel
      
      
      # sidebarPanel(
      #   tags$h3('graph settings'),
      #   sliderInput('rotation', 'rotate labels', 0, 5, 0, 1)
      # ),
      # 
      
      mainPanel(
        #h1("Header 1"),
        
        #h4("Output 1"),
        #verbatimTextOutput("txtout"),
        plotOutput('lm_plot', brush = "plot_brush",width = '4000px'),
        tableOutput("lm_plot_info"),
        plotOutput('luca_plot', brush = "plot_brush", width = 'auto'),
        tableOutput("luca_plot_info"),
        
        
        #h4('Read table'),
        tableOutput('edds'),
        tableOutput('flowjo'),
        tableOutput('dosing'),
        tableOutput('combined'),
        
      ) # mainPanel
      
    ),
    # Navbar 1, tabPanel

        
    tabPanel("TMDD", 
             sidebarPanel(
               tags$h3("Convert MFI to ABC"),
               fileInput("quickcal", "Upload filled quick cal .xls file", accept = c(".xls",'.xlsx')),
               fileInput('quickcal_conv_list', 'Upload .csv file(s) to be converted into ABC', accept = c('.csv'), multiple = TRUE),
               downloadButton("download_abc", "Download converted files"),
               
             
               tags$h3("Calculations using EDDS file"),
               fileInput("edds_abc", "Upload EDDS .csv", accept = ".csv"),
               fileInput("bead_fls_mabs", "Upload .csv file containing mfi of beads with corresponding mAb as header", accept = ".csv"),
               fileInput(
                 "flowjo_abc",
                 "Upload flowjo saved .xlsx files",
                 accept = c('.xls', ".xlsx"),
                 multiple = TRUE
               ),
               textInput('mfi_var_abc','Enter MFI variable', 'Geometric Mean : pHAb-A') , textInput('qe','Enter calculated %QE')
             ),
             
             mainPanel(
               tableOutput('test_abc')
               
               
             )
             
             
             ),
  
  
  
    tabPanel("DCIA", "This panel is intentionally left blank"),
    tabPanel("SEC", "This panel is intentionally left blank"),
    # tabPanel("EDDS_STITCH", "This panel is intentionally left blank"),
    tabPanel("MIXING"),
    tabPanel("WORKFLOW", "This panel is intentionally left blank",
    
             imageOutput('workflow', width = 'auto')
             
             
             
             
             ),
    tabPanel("VISUALIZE",
             sidebarLayout(
               sidebarPanel(
                 radioButtons(
                   'vizoptions',
                   label = 'Select Plate',
                   choices = c('plate1', 'plate2'),
                   selected = 'plate1'
                 )
               ),
               
               mainPanel(plotOutput('viz'))
               
               
             )),
    
    
    tabPanel("HOMER",
             
             sidebarLayout(
               sidebarPanel(
                 imageOutput('homer', width = 'auto'),
                 radioButtons(
                   'lang',
                   label = NULL,
                   choices = c('English', 'German'),
                   selected = 'English'
                 ),
                 actionButton('randomise', 'Randomise')
                 
                 
                 
               ),
               
               mainPanel(htmlOutput("quote"),) # mainPanel
               
             ))
    
  ) # navbarPage
) # fluidPage


# Define server function
server <- function(input, output, session) {
  
  edds <- reactive({
    req(input$edds)
    
    edds_process(input$edds$datapath)
  })
  
  
  # flowjo file processing --------------------------------------------------
  
  
  
  flowjo <- reactive({
    req(input$flowjo)
    flowjo_processing(input$flowjo$datapath)
  })
  
  
  
  # observe({
  #   updateSelectInput(session, "mfi_choices",
  #                     choices = colnames(flowjo())
  # )})
# Tecan fluorescence dosing processing ------------------------------------
  
  
  
  dosing <- reactive({
    req(input$dosing)
    dosing_processing(input$dosing$datapath)
    })
  
  
  # combine EDDS with dosing and flowjo data --------------------------------
  
  
  mfi_choices <- reactive(input$mfi_choices)
  
  control_mabs <- reactive(input$control_mabs |> 
                             str_split(',', simplify = TRUE) |>
                             str_trim())
  
  
  
  EDDS_combined <- reactive({
    
    req(edds(), mfi_choices(), flowjo(), dosing())
    EDDS_combined_processing(edds(),mfi_choices(),flowjo(),dosing())
    
    
  })
      
  
  # Mathematical manipulation on EDDS() -------------------------------------

  EDDS_combined_analysed <- reactive({
    req(EDDS_combined())
    req(mfi_choices())
    req(control_mabs())
    
    edds_dn <- EDDS_combined()
    
    edds_dn <-  edds_dn |> 
      mutate(`Viability` = `Cell count_morphology_live` / `Cell count_morphology`,
             `Cell fraction_gated` = `Cell count_morphology_live` /`Cell count_total`) |> 
      group_by(`Experiment date`, `Biosample ID`, `Incubation time`) %>%
      do(mock_sub(.,mfi_choices())) |>
      mutate(
        '{mfi_choices()}_BG subtracted_dose normalized' := .data[[paste0(mfi_choices(), '_BG subtracted')]] /
          `Fluorescence_dosing solution`
      ) 
    
    
    edds_dn %>% filter(!str_detect(`Tapir ID_unlabeled molecule (parent)`,
                                   regex('mock', ignore_case = TRUE))) %>%
      group_by(`Experiment date`,
               `Tapir ID_unlabeled molecule (parent)`,
               `Biosample ID`) %>%
      do(lm_wo_normpoint = lm(.[[paste0(mfi_choices(), '_BG subtracted_dose normalized')]] ~ 0 + .[['Incubation time']],data = .)
      ) -> model_wo.pointnorm
    
    edds_dn <- open_model(model_wo.pointnorm,'lm_wo_normpoint',edds_dn)
    
    
    edds_dn <- edds_dn %>% ungroup() %>% group_by(`Experiment date`,`Biosample ID`) %>%
      do(min_max_day(.,control_mabs(),'estimate__lm_wo_normpoint')) %>% ungroup()
    
    
    
    edds_dn <-
      edds_dn %>% 
      group_by(`Experiment date`, `Biosample ID`) %>% 
      do(point_norm(., mfi_choices(), control_mabs())) %>% ungroup()
    
    model_w.pointnorm <- edds_dn %>% filter(!str_detect(`Tapir ID_unlabeled molecule (parent)`,
                                                        regex('mock', ignore_case = TRUE))) %>% 
      group_by(`Tapir ID_unlabeled molecule (parent)`,`Biosample ID`) %>%
      do(lm_w_normpoint = lm(.[[paste0(mfi_choices(), '_BG subtracted_dose_control normalized')]] ~ 0 + .[['Incubation time']],data = .))
    
    edds_dn <- open_model(model_w.pointnorm,'lm_w_normpoint',edds_dn)
    

  })
  
  
  
  
  
  
  
  
  # homer -------------------------------------------------------------------
  
  output$homer <- renderImage({
    list(src = r'(homer.jpg)')
  }, deleteFile = FALSE)
  
  lang_path <- reactive(switch(input$lang,
                               'English' = r'(homer_quotes.txt)',
                               'German' = r'(homer_quotes_german.txt)'))
  
  
  homer_quote <- eventReactive(input$randomise, {
    sample(read_lines(lang_path()),
           size = 1)
  })
  
  
  
  output$quote <- renderUI({
    # quote_english <- read_lines(r'(C:\Users\PADAMSEA\Desktop\r app scripts\homer\homer_quotes.txt)')
    HTML(paste0("<h2>", homer_quote(), "</h2>"))
  })
  
  
  
  # plate viz ---------------------------------------------------------------
  
  
  
  output$viz <- renderPlot({
    edds() |> group_split(`Plate number`) %>% .[[1]] |> select(`Well number`, `Tapir ID_unlabeled molecule (parent)`) |>
      
      separate(
        `Well number`,
        sep = 1,
        into = c('pl_row', 'pl_col'),
        convert = TRUE
      ) -> plate_layout
    
    plate_layout_df <-
      tibble(pl_col = rep(x, each = 8), pl_row = rep(y, 12))
    plot_join <- full_join(plate_layout, plate_layout_df)
    
    # test_join |> mutate(V2 = case_when(V2,is.na(),'PBS'))
    
    ggplot(plot_join) +
      geom_point(size = 8,
                 aes(pl_col, pl_row, color = `Tapir ID_unlabeled molecule (parent)`)) +
      theme_classic() +
      scale_y_discrete(limits = rev) +  #scale_y_discrete(labels= rev(y))
      scale_x_continuous(breaks = x) +
      xlab('') +
      ylab('')
    
    
  })
  
  
  
  # PREVIEWS ----------------------------------------------------------------
  
  

  output$edds <- renderTable({
     head(edds())

  })

  output$flowjo <- renderTable({
    head(flowjo())
  })

  output$dosing <- renderTable({
    head(dosing())
  })

  output$combined <- renderTable({
    head(EDDS_combined())
  })



# plots -------------------------------------------------------------------
# 
  output$lm_plot <- renderPlot({
    req(EDDS_combined_analysed())

    EDDS_combined_analysed() |> ggplot(aes(`Incubation time`,`Geometric Mean : pHAb-A_BG subtracted_dose normalized`)) +
      geom_smooth(method='lm', formula = 'y ~ 0 + x', fullrange=TRUE) +
      geom_point() +
      # facet_wrap(vars(`Tapir ID_unlabeled molecule (parent)`)) + 
      expand_limits(x = 0, y = 0) +
      facet_grid(rows=vars(`Experiment date`),cols = vars(`Tapir ID_unlabeled molecule (parent)`))

  })
#   
#   
#   output$lm_plot_info <- renderTable({
#     req(EDDS_combined_analysed())
#     brushedPoints(test_edds, input$plot_brush)
#   })
#   
#   
#   
#   
#   output$luca_plot <- renderPlot({
#     req(EDDS_combined_analysed())
#     test_edds |> ggplot(aes(`Tapir ID_unlabeled molecule (parent)`,slope)) +
#       geom_boxplot() + geom_jitter() 
#     
#   })
#   
#   output$luca_plot_info <- renderTable({
#     req(EDDS_combined_analysed())
#     brushedPoints(test_edds, input$plot_brush)
#   })
#     
    
  
  
  
  
  # Misc --------------------------------------------------------------------
  # output$misc <- reactive(mfi_choices())
  
  
  output$control_choices <- reactive({
    req(input$edds)
    edds()$`Tapir ID_unlabeled molecule (parent)` |> unique() |> c() |> paste0()
    
  })
  
  output$mfi_choices_list <- reactive({
    req(flowjo())
    names(flowjo())[(which(flowjo() |> names()  == 'ultimate_gate') + 1):length(flowjo())]
  })
  
  
  
  output$workflow <- renderImage({
    list(src = 'rshiny.png')
  }, deleteFile = FALSE)
  
  
  

# TMDD --------------------------------------------------------------------
  
 quickcal <- reactive({
  req(input$quickcal)
   quickcal <-  readxl::read_excel(input$quickcal$datapath)[7:10,3:4] %>% mutate(across(everything(),as.numeric))
   colnames(quickcal) <- c('abc','bead_fl')  
   return(quickcal)                                             
                                                })
  
  
 to_abc <- reactive({
   req(input$quickcal_conv_list)
   
   map(input$quickcal_conv_list$datapath,\(x) read_csv(x,col_names = FALSE) %>% 
         mutate(across(everything(),as.numeric)))
   
  })
 
 
   
 converted_abc_list <- reactive({
   
   req(to_abc(),quickcal())
   
   model_abc <- lm(log(abc) ~ log(bead_fl), data = quickcal()) #make linear model (on log scale)
   
   to_abc() %>% map(
     
     \(x) x %>% map(
       
       \(y) predict(model_abc,newdata = list('bead_fl' = y)) %>% exp()
       
     )
   ) %>% map(as.data.frame)
   
   
 })
 
 
 
   
 output$test_abc <- renderTable(converted_abc_list()[[1]])
  
 output$download_abc <- downloadHandler(
   filename = function(){
     paste("converted_abc", Sys.Date(), ".zip", sep = "")
   },
   content = function(file){
     
     temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
     dir.create(temp_directory)
     
     converted_abc_list() %>%
       map2(tools::file_path_sans_ext(input$quickcal_conv_list$name),function(x,y){
         if(!is.null(x)){
           file_name <- glue("{y}_abc.csv")
           readr::write_csv(x, file.path(temp_directory, file_name),col_names = FALSE)
         }
       })
     
     
     zip::zip(
       zipfile = file,
       files = dir(temp_directory),
       root = temp_directory
     )
     
     
     
   },
   contentType = "application/zip"
   
 )
  
  
  
  
# DOWNLOAD ----------------------------------------------------------------
  
  output$download_pzfx <-  downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$edds$name), ".pzfx")
    },
    content = function(file) {
      pzfx::write_pzfx(
        graphpad(EDDS_combined_analysed(),mfi_choices()), file, x_col = 1)
      
    })
  
  
 
  output$download_csv <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$edds$name), ".csv")
    },
    content = function(file) {
      write_csv(EDDS_combined_analysed(),file)
  })
  # output$contents <- renderTable({
  # #   file <- input$two_hr
  # #   ext <- tools::file_ext(file$datapath)
  # #
  #   # req(file)
  #   # validate(need(ext == ".xls", "Please upload a xls or xlsx file"))
  #   #
  #   two_hr_read <- readxl::read_excel(file$datapath)
  #   #write.csv(two_hr_read,'C:\\Users\\PADAMSEA\\Downloads\\test.csv')
  #   clean_list(two_hr_read)[[1]]
  # })
} # server


# Create Shiny object
shinyApp(ui = ui, server = server)







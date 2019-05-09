library(shiny)
library(taxize)

# Define UI for app that draws a histogram ----

kingdom_list <- c("Archaea","Bacteria","Eukaryota","Fungi","Metazoa","Plant"="Viridiplantae","Viruses")

ui <- fluidPage(
  
  fluid=TRUE,
  theme = "bootstrap.min.css",
  
  # App title ----
  titlePanel("Hello Shiny!"),
  
  tabPanel(title = "Library Generation",
           mainPanel(# the following lines could be uncommented when the download ref seq can
             # work on the rest of the kingdoms
             # radioButtons("kingdom", "Choose a kingdom:",
             #              c("Archaea" = "archaea",
             #                "Bacteria" = "bacteria",
             #                "Fungi" = "fungi",
             #                "Invertebrate" = "invertebrate",
             #                "Plant" = "plant",
             #                "Protozoa" = "protozoa",
             #                "Vertebrate" = "vertibrate",
             #                "Vertebrate other" = "vertibrate_other",
             #                "Virus" = "viral")
             #              ),
             
             # create checkbox input for representative library and reference library
             checkboxInput("representative", "representative", value = TRUE, width = NULL),
             checkboxInput("reference", "reference", value = FALSE, width = NULL),
             
             
             actionButton("downloadref","Download Ref_Seq")
           )
  ),
  mainPanel(
    tabsetPanel(
      type = "tabs",
      
      tabPanel("Kingdom",
               radioButtons("kingdomGroup", label = "Choose a kingdom:", 
                                  choices = kingdom_list),
               actionButton("kingdom_update","Update")
      ),
      
      tabPanel("Phylum",
               checkboxGroupInput("phylumGroup", label = "Choose a phylum:"),
               actionLink("phylum_selectall","select all"),
               actionButton("phylum_update","Update")
      ),
      
      tabPanel("Class",
               checkboxGroupInput("classGroup", label = "Choose a class:"),
               actionLink("class_selectall","select all"),
               actionButton("class_update","Update")
      ),
      tabPanel("Order",
               checkboxGroupInput("orderGroup", label = "Choose an order:"),
               actionLink("order_selectall","select all"),
               actionButton("order_update","Update")
      ),
      tabPanel("Family",
               checkboxGroupInput("familyGroup", label = "Choose a family:"),
               actionLink("family_selectall","select all"),
               actionButton("family_update","Update")
      ),
      tabPanel("Genus",
               checkboxGroupInput("genusGroup", label = "Choose a genus:"),
               actionLink("genus_selectall","select all"),
               actionButton("genus_update","Update")
      ),
      tabPanel("Species",
               checkboxGroupInput("speciesGroup", label = "Choose a species:"),
               actionLink("species_selectall","select all"),
               actionButton("species_update","Update")
      )
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  output$value <- renderPrint({ input$kingdomGroup })
  
  #output$phylumPanel <- renderUI({
  #  if(!is.null(input$kingdomGroup)){
  #    checkboxGroupInput("phylumGroup", label = "Choose a phylum:", 
  #                       choices = c("A","B"))
  #  }else{
  #    h3("No kingdom selected")
  #  }
   # 
  #})
  
  p_list.k <- reactiveVal(character(0))
  c_list.k <- reactiveVal(character(0))
  o_list.k <- reactiveVal(character(0))
  f_list.k <- reactiveVal(character(0))
  g_list.k <- reactiveVal(character(0))
  s_list.k <- reactiveVal(character(0))
  
  c_list.p <- reactiveVal(character(0))
  o_list.p <- reactiveVal(character(0))
  f_list.p <- reactiveVal(character(0))
  g_list.p <- reactiveVal(character(0))
  s_list.p <- reactiveVal(character(0))
  
  o_list.c <- reactiveVal(character(0))
  f_list.c <- reactiveVal(character(0))
  g_list.c <- reactiveVal(character(0))
  s_list.c <- reactiveVal(character(0))
  
  f_list.o <- reactiveVal(character(0))
  g_list.o <- reactiveVal(character(0))
  s_list.o <- reactiveVal(character(0))
  
  g_list.f <- reactiveVal(character(0))
  s_list.f <- reactiveVal(character(0))
  
  s_list.g <- reactiveVal(character(0))
  
  
  observeEvent(input$kingdom_update,{
    k_input <- input$kingdomGroup
    
    p_list.new <- character(0)
    c_list.new <- character(0)
    o_list.new <- character(0)
    f_list.new <- character(0)
    g_list.new <- character(0)
    s_list.new <- character(0)
    
    for(k in k_input){
      children_list <- NULL
      children_list <- children(k, db = 'ncbi')
      if(!is.null(children_list)){
        children_list <- children_list[[1]]
        
        children_list.p <- children_list[which(children_list$childtaxa_rank == 'phylum'),]
        children_list.c <- children_list[which(children_list$childtaxa_rank == 'class'),]
        children_list.o <- children_list[which(children_list$childtaxa_rank == 'order'),]
        children_list.f <- children_list[which(children_list$childtaxa_rank == 'family'),]
        children_list.g <- children_list[which(children_list$childtaxa_rank == 'genus'),]
        children_list.s <- children_list[which(children_list$childtaxa_rank == 'species'),]
        
        p_list.new <- append(p_list.new, children_list.p$childtaxa_name)
        c_list.new <- append(c_list.new, children_list.c$childtaxa_name)
        o_list.new <- append(o_list.new, children_list.o$childtaxa_name)
        f_list.new <- append(f_list.new, children_list.f$childtaxa_name)
        g_list.new <- append(g_list.new, children_list.g$childtaxa_name)
        s_list.new <- append(s_list.new, children_list.s$childtaxa_name)
      }
    }
    
    updateCheckboxGroupInput(session, "phylumGroup", choices = sort(p_list.new))
    updateCheckboxGroupInput(session, "classGroup", choices = sort(c_list.new))
    updateCheckboxGroupInput(session, "orderGroup", choices = sort(o_list.new))
    updateCheckboxGroupInput(session, "familyGroup", choices = sort(f_list.new))
    updateCheckboxGroupInput(session, "genusGroup", choices = sort(g_list.new))
    updateCheckboxGroupInput(session, "speciesGroup", choices = sort(s_list.new))
    
    p_list.k(p_list.new)
    c_list.k(c_list.new)
    o_list.k(o_list.new)
    f_list.k(f_list.new)
    g_list.k(g_list.new)
    s_list.k(s_list.new)
  })
  
  observeEvent(input$phylum_update,{
    p_input <- input$phylumGroup
    
    c_list.new <- character(0)
    o_list.new <- character(0)
    f_list.new <- character(0)
    g_list.new <- character(0)
    s_list.new <- character(0)
    
    withProgress(message="Getting children taxons from selected phyla",value=0,{
      for(p in p_input){
        children_list <- NULL
        Sys.sleep(0.1)
        children_list <- children(p, db = 'ncbi')
        if(!is.null(children_list)){
          children_list <- children_list[[1]]
          
          children_list.c <- children_list[which(children_list$childtaxa_rank == 'class'),]
          children_list.o <- children_list[which(children_list$childtaxa_rank == 'order'),]
          children_list.f <- children_list[which(children_list$childtaxa_rank == 'family'),]
          children_list.g <- children_list[which(children_list$childtaxa_rank == 'genus'),]
          children_list.s <- children_list[which(children_list$childtaxa_rank == 'species'),]
          
          c_list.new <- append(c_list.new, children_list.c$childtaxa_name)
          o_list.new <- append(o_list.new, children_list.o$childtaxa_name)
          f_list.new <- append(f_list.new, children_list.f$childtaxa_name)
          g_list.new <- append(g_list.new, children_list.g$childtaxa_name)
          s_list.new <- append(s_list.new, children_list.s$childtaxa_name)
        }
        incProgress(1/length(p_input))
      }
    })
    
    updateCheckboxGroupInput(session, "classGroup", choices = sort(unique(c(c_list.new,
                                                                            c_list.k()
                                                                            ))))
    updateCheckboxGroupInput(session, "orderGroup", choices = sort(unique(c(o_list.new,
                                                                            o_list.k()
                                                                            ))))
    updateCheckboxGroupInput(session, "familyGroup", choices = sort(unique(c(f_list.new,
                                                                             f_list.k()
                                                                             ))))
    updateCheckboxGroupInput(session, "genusGroup", choices = sort(unique(c(g_list.new,
                                                                            g_list.k()
                                                                            ))))
    updateCheckboxGroupInput(session, "speciesGroup", choices = sort(unique(c(s_list.new,
                                                                              s_list.k()
                                                                              ))))
    
    c_list.p(c_list.new)
    o_list.p(o_list.new)
    f_list.p(f_list.new)
    g_list.p(g_list.new)
    s_list.p(s_list.new)
  })
  
  observeEvent(input$class_update,{
    c_input <- input$classGroup
    
    o_list.new <- character(0)
    f_list.new <- character(0)
    g_list.new <- character(0)
    s_list.new <- character(0)
    
    withProgress(message="Getting children taxons from selected classes", value=0, {
        for(c in c_input){
          children_list <- NULL
          Sys.sleep(0.1)
          children_list <- children(c, db = 'ncbi')
          if(!is.null(children_list)){
            children_list <- children_list[[1]]
            
            children_list.o <- children_list[which(children_list$childtaxa_rank == 'order'),]
            children_list.f <- children_list[which(children_list$childtaxa_rank == 'family'),]
            children_list.g <- children_list[which(children_list$childtaxa_rank == 'genus'),]
            children_list.s <- children_list[which(children_list$childtaxa_rank == 'species'),]
            
            o_list.new <- append(o_list.new, children_list.o$childtaxa_name)
            f_list.new <- append(f_list.new, children_list.f$childtaxa_name)
            g_list.new <- append(g_list.new, children_list.g$childtaxa_name)
            s_list.new <- append(s_list.new, children_list.s$childtaxa_name)
          }
          incProgress(1/length(c_input))
        }
    })
    
    updateCheckboxGroupInput(session, "orderGroup", choices = sort(unique(c(o_list.new,
                                                                            o_list.k(),
                                                                            o_list.p()
                                                                            ))))
    updateCheckboxGroupInput(session, "familyGroup", choices = sort(unique(c(f_list.new,
                                                                             f_list.k(),
                                                                             f_list.p()
                                                                             ))))
    updateCheckboxGroupInput(session, "genusGroup", choices = sort(unique(c(g_list.new,
                                                                            g_list.k(),
                                                                            g_list.p()
                                                                            ))))
    updateCheckboxGroupInput(session, "speciesGroup", choices = sort(unique(c(s_list.new,
                                                                              s_list.k(),
                                                                              s_list.p()
                                                                              ))))
    
    o_list.c(o_list.new)
    f_list.c(f_list.new)
    g_list.c(g_list.new)
    s_list.c(s_list.new)
  })
  
  observeEvent(input$order_update,{
    o_input <- input$orderGroup

    f_list.new <- character(0)
    g_list.new <- character(0)
    s_list.new <- character(0)
    
    withProgress(message='Getting children taxons from selected orders',value=0,{
      for(o in o_input){
        children_list <- NULL
        Sys.sleep(0.1)
        children_list <- children(o, db = 'ncbi')
        if(!is.null(children_list)){
          children_list <- children_list[[1]]
          
          children_list.f <- children_list[which(children_list$childtaxa_rank == 'family'),]
          children_list.g <- children_list[which(children_list$childtaxa_rank == 'genus'),]
          children_list.s <- children_list[which(children_list$childtaxa_rank == 'species'),]
          
          f_list.new <- append(f_list.new, children_list.f$childtaxa_name)
          g_list.new <- append(g_list.new, children_list.g$childtaxa_name)
          s_list.new <- append(s_list.new, children_list.s$childtaxa_name)
        }
        incProgress(1/(length(o_input)))
      }
    })
    
    updateCheckboxGroupInput(session, "familyGroup", choices = sort(unique(c(f_list.new,
                                                                             f_list.k(),
                                                                             f_list.p(),
                                                                             f_list.c()
                                                                             ))))
    updateCheckboxGroupInput(session, "genusGroup", choices = sort(unique(c(g_list.new,
                                                                            g_list.k(),
                                                                            g_list.p(),
                                                                            g_list.c()
                                                                            ))))
    updateCheckboxGroupInput(session, "speciesGroup", choices = sort(unique(c(s_list.new,
                                                                              s_list.k(),
                                                                              s_list.p(),
                                                                              s_list.c()
                                                                              ))))
    
    f_list.o(f_list.new)
    g_list.o(g_list.new)
    s_list.o(s_list.new)
  })
  
  observeEvent(input$family_update,{
    f_input <- input$familyGroup
    
    g_list.new <- character(0)
    s_list.new <- character(0)
    
    withProgress(message='Getting children taxon from selected families', value=0,{
      for(f in f_input){
        children_list <- NULL
        Sys.sleep(0.1)
        children_list <- children(f, db = 'ncbi')
        if(!is.null(children_list)){
          children_list <- children_list[[1]]
          
          children_list.g <- children_list[which(children_list$childtaxa_rank == 'genus'),]
          children_list.s <- children_list[which(children_list$childtaxa_rank == 'species'),]
  
          g_list.new <- append(g_list.new, children_list.g$childtaxa_name)
          s_list.new <- append(s_list.new, children_list.s$childtaxa_name)
        }
        incProgress(1/length(f_input))
      }
    })
    
    updateCheckboxGroupInput(session, "genusGroup", choices = sort(unique(c(g_list.new,
                                                                            g_list.k(),
                                                                            g_list.p(),
                                                                            g_list.c(),
                                                                            g_list.o()
                                                                            ))))
    updateCheckboxGroupInput(session, "speciesGroup", choices = sort(unique(c(s_list.new,
                                                                              s_list.k(),
                                                                              s_list.p(),
                                                                              s_list.c(),
                                                                              s_list.o()
                                                                              ))))
    
    g_list.f(g_list.new)
    s_list.f(s_list.new)
  })
  
  observeEvent(input$genus_update,{
    g_input <- input$genusGroup
    
    s_list.new <- character(0)
    
    withProgress(message='Getting children taxons from selected genuses',value=0,{
      for(g in g_input){
        children_list <- NULL
        Sys.sleep(0.5)
        children_list <- children(g, db = 'ncbi')
        if(!is.null(children_list)){
          children_list <- children_list[[1]]
          
          children_list.s <- children_list[which(children_list$childtaxa_rank == 'species'),]
          
          s_list.new <- append(s_list.new, children_list.s$childtaxa_name)
        }
        incProgress(1/length(g_input))
      }
    })
    
    updateCheckboxGroupInput(session, "speciesGroup", choices = sort(unique(c(s_list.new,
                                                                              s_list.k(),
                                                                              s_list.p(),
                                                                              s_list.c(),
                                                                              s_list.o(),
                                                                              s_list.f()
                                                                              ))))
    
    s_list.g(s_list.new)
  })
  
  # set action for "select all" actionLinks
  observeEvent(input$phylum_selectall,{
    if(input$phylum_selectall%%2 == 0){
      updateCheckboxGroupInput(session,
                               "phylumGroup",
                               choices=sort(unique(p_list.k()))
      )
    }
    else{
      updateCheckboxGroupInput(session,
                               "phylumGroup",
                               choices=sort(unique(p_list.k())),
                               selected=sort(unique(p_list.k()))
      )
    }
  })
  observeEvent(input$class_selectall,{
    if(input$class_selectall%%2 == 0){
      updateCheckboxGroupInput(session,
                               "classGroup",
                               choices=sort(unique(c(c_list.k(),c_list.p())))
      )
    }
    else{
      updateCheckboxGroupInput(session,
                               "classGroup",
                               choices=sort(unique(c(c_list.k(),c_list.p()))),
                               selected=sort(unique(c(c_list.k(),c_list.p())))
      )
    }
  })
  observeEvent(input$order_selectall,{
    if(input$order_selectall%%2 == 0){
      updateCheckboxGroupInput(session,
                               "orderGroup",
                               choices=sort(unique(c(o_list.k(),o_list.p(),o_list.c())))
      )
    }
    else{
      updateCheckboxGroupInput(session,
                               "orderGroup",
                               choices=sort(unique(c(o_list.k(),o_list.p(),o_list.c()))),
                               selected=sort(unique(c(o_list.k(),o_list.p(),o_list.c())))
      )
    }
  })
  observeEvent(input$family_selectall,{
    if(input$family_selectall%%2 == 0){
      updateCheckboxGroupInput(session,
                               "familyGroup",
                               choices=sort(unique(c(f_list.k(),f_list.p(),f_list.c(),f_list.o())))
      )
    }
    else{
      updateCheckboxGroupInput(session,
                               "familyGroup",
                               choices=sort(unique(c(f_list.k(),f_list.p(),f_list.c(),f_list.o()))),
                               selected=sort(unique(c(f_list.k(),f_list.p(),f_list.c(),f_list.o())))
      )
    }
  })
  observeEvent(input$genus_selectall,{
    if(input$genus_selectall%%2 == 0){
      updateCheckboxGroupInput(session,
                               "genusGroup",
                               choices=sort(unique(c(g_list.k(),g_list.p(),g_list.c(),g_list.o(),g_list.f())))
      )
    }
    else{
      updateCheckboxGroupInput(session,
                               "genusGroup",
                               choices=sort(unique(c(g_list.k(),g_list.p(),g_list.c(),g_list.o(),g_list.f()))),
                               selected=sort(unique(c(g_list.k(),g_list.p(),g_list.c(),g_list.o(),g_list.f())))
      )
    }
  })
  observeEvent(input$species_selectall,{
    if(input$species_selectall%%2 == 0){
      updateCheckboxGroupInput(session,
                               "speciesGroup",
                               choices=sort(unique(c(s_list.k(),s_list.p(),s_list.c(),s_list.o(),s_list.f(),s_list.g())))
      )
    }
    else{
      updateCheckboxGroupInput(session,
                               "speciesGroup",
                               choices=sort(unique(c(s_list.k(),s_list.p(),s_list.c(),s_list.o(),s_list.f(),s_list.g()))),
                               selected=sort(unique(c(s_list.k(),s_list.p(),s_list.c(),s_list.o(),s_list.f(),s_list.g())))
      )
    }
  })
  
  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)

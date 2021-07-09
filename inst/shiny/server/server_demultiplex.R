#Barcode File Display

#Inputting data
reactbarcode <- reactive ({
  req(input$barcode)
  barcodeFile <-
    read.table(input$barcode$datapath,
               header = input$headerbarcode,
               sep = "\t")
  barcodeFile
})
#Updates slider
observe({
  updateSliderInput (
    session,
    "displaybarcode",
    min = 1,
    value = 10,
    max = nrow(reactbarcode())
  )
})
#Display Index Data
updatedisplaybarcode <- reactive({
  displayb <- reactbarcode()[1:input$displaybarcode, ]
  return(displayb)
})
#Output Data
output$tablebarcode <- renderDataTable({
  updatedisplaybarcode()[1:input$displaybarcode, ]
})

#Update Hamming Distance Maximum
observe({
  updateSelectInput(session,
                    "hammingDist",
                    choices = c(1:nchar(as.character(reactbarcode(
                    )[1, 2]))),
                    selected = 1)
})

#Index and Read File Display

#Inputting Index data
reactindex <- reactive ({
  req(input$index)
  indexFile <-
    Biostrings::readDNAStringSet(input$index$datapath, format = "fastq", n =
                                   100)
  indexFile
})

#Inputting Read data
reactread <- reactive ({
  req(input$read)
  readFile <-
    Biostrings::readQualityScaledDNAStringSet(input$read$datapath, n = 100)
  readFile
})

#Updates slider
observe({
  updateSliderInput (session,
                     "display",
                     min = 1,
                     value = 10,
                     max = 100)
})
observeEvent(input$upload, {
  #Display Data
  updatedisplay <- reactive({
    Index <- reactindex()[1:input$display, ]
    Read <- reactread()[1:input$display, ]
    x <- cbind.data.frame(Index, Read)
    return(x)
  })
  #Output Data
  output$table <- renderTable({
    updatedisplay()[1:input$display, ]
  })
})


#rcBarcodes

rc <- reactive({
  return(as.logical(input$rcBarcodes))
})

#hammingDist

ham <- reactive({
  return(as.numeric(input$hammingDist))
})


out <- reactive({
  return(as.character(input$location))
})

#Demultiplex
observeEvent(input$demultiplex, {
  demultiplex(
    barcodeFile = input$barcode$datapath ,
    indexFile = input$index$datapath ,
    readFile = input$read$datapath,
    location = file.path(out(), "/Demultiplexed Sample Files"),
    rcBarcodes = rc(),
    cores = cores(),
    hammingDist = ham()
  )
  
})





  
  
  
  
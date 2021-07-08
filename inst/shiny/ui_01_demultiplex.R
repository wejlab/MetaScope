
extractReads <- function(barcodeIndex, barcodes, sampleNames, index, reads, 
                         location = "./demultiplex_fastq", rcBarcodes = TRUE, hDist = 0) {
  barcode <- barcodes[barcodeIndex]
  sampleName <- sampleNames[barcodeIndex]
  message("Finding reads for barcode: ", barcode)
  if (rcBarcodes) {
    rci <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(barcode)))
  } else {
    rci <- barcode
  }
  # ind_match <- as.character(index) == rci
  ind_match <- adist(as.character(index), rci) <= hDist
  
  numReads <- sum(ind_match)
  outFileName <- paste(location, "/", sampleName, "_", barcode, ".fastq.gz", 
                       sep = "")
  if (numReads == 0) {
    message("\tFound 0 reads for this barcode, no file will be written")
  } else {
    message("\tFound ", sum(ind_match), " reads, writing reads to: ", 
            outFileName)
    Biostrings::writeQualityScaledXStringSet(reads[c(ind_match)], outFileName, 
                                             compress = TRUE)
  }
  return(list(output_file = outFileName, numberOfReads = numReads, matchedIndexes = ind_match))
}

    
demultiplex <- function(bcFile, inds, reads, rcBarcodes=TRUE, 
                        location = file.path(location,"/demultiplex_fastq"), cores=1,    hammingDist=0) {
  message("Reading Sample Names and Barcodes from: ", bcFile)
  barcodes <- bcFile[, 2]
  samNames <- bcFile[, 1]
  message("\tFound information for ", length(barcodes), " samples/barcodes")
  
  message("Reading Index File: ", inds)
  message("\tFound indexes for ", length(inds), " reads")
  
  message("Reading Sequence File: ", reads)
  message("\tFound ", length(reads), " reads")
  
  ## make output directory if nessary
  if (!dir.exists(location)) {
    dir.create(location)
  }
  
  # Loop over barcodes
  numReads <- NULL
  ind_no_match <- numeric(length(reads))
  for (i in 1:length(barcodes)) {
    extracted <- extractReads(i, barcodes, samNames, inds, reads, rcBarcodes = rcBarcodes, 
                              location = location, hDist = hammingDist)
    numReads <- c(numReads, extracted$numberOfReads)
    ind_no_match <- ind_no_match + extracted$matchedIndexes
  }
  message(sum(ind_no_match > 1))
  ind_no_match <- (ind_no_match == 0)
  
  # sapply over barcodes and writing to file -- for multi-threading if
  # (cores == 1){ extracted <- lapply(1:length(barcodes), extractReads,
  # barcodes, samNames, inds, reads, hammingDist) }else{ message('Using
  # ',cores,' cores') multicoreParam <- MulticoreParam(workers = cores)
  # extracted <- bplapply(1:length(barcodes),extractReads, barcodes,
  # samNames, inds, reads, BPPARAM = multicoreParam) } numReads <-
  # sapply(extracted, function(x) x$numberOfReads) ind_no_match <- (
  # rowSums(sapply(extracted, function(x) x$matchedIndexes)) == 0 )
  
  # number of reads for each barcode
  if (any(numReads == 0)) {
    message("Did not find any reads for the following barcodes: ", 
            paste(barcodes[numReads == 0], collapse = " "))
    message("Did not find any reads for the following samples: ", paste(samNames[numReads == 
                                                                                   0], collapse = " "))
    write(paste("Did not find any reads for the following barcodes:", 
                paste(barcodes[numReads == 0], collapse = " "), "\n", "Did not find any reads for the following samples: ", 
                paste(samNames[numReads == 0], collapse = " ")), file = "/demultiplex_fastq/unmapped_barcodes_samples.txt")
  }
  
  # Track reads without matches, and write them to an 'orphan' file
  message(paste("Found ", sum(ind_no_match), " reads without a matching barcode (", 
                100 * round(mean(ind_no_match), 4), "%), writing reads to: ", location, 
                "/orpahns.fastq.gz", sep = ""))
  Biostrings::writeQualityScaledXStringSet(reads[c(ind_no_match)], paste(location, 
                                                                         "/orpahns.fastq.gz", sep = ""), compress = TRUE)
  
  summaryMat <- cbind(bcFile[1:length(barcodes), ], NumberOfReads = numReads)
  write.table(summaryMat, file = paste(location, "/summary.txt", sep = ""), 
              col.names = FALSE, row.names = TRUE, quote = FALSE)
  return(summaryMat)
} 






library(shiny)
library(DT)
library(Biostrings)
library(rsconnect)

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
#Read in Barcode File 
         
         fileInput("barcode", "Barcode File",
                accept = c(
                  ".txt")),
         checkboxInput("headerbarcode", "Header", TRUE),
         sliderInput("displaybarcode","View Rows:",min=0,max=0,value=10,step=1),
         
#Read in Read File          

         fileInput("read", "Read File",
                accept = c(
                   ".fast1.gz")),
         sliderInput("displayread","View Rows:",min=0,max=0,value=10,step=1),

#Read in Index File
      
         fileInput("index", "Index File",
                accept = c(
                  ".fast1.gz")),
         sliderInput("displayindex","View Rows:",min=0,max=0,value=10,step=1),

#Hamming, rcBarcode, and Codes

      selectInput("hammingDist","Hamming",choices=c(1:20),selected=1),
      checkboxInput("rcBarcodes", "rcBarcode", TRUE),
      selectInput("cores","Cores",choices=c(1:10),selected=1),

#Specify Output File

       textInput("location", "Output File Location",
                  placeholder = "New Folder"
      ),

#Dimultiplex!

      actionButton("demultiplex","Demultiplex")
    ),

#Index Table Display 

    mainPanel(
      tableOutput("tableindex"),tableOutput("tableread"), dataTableOutput("tablebarcode")
    )
   )
  )

server <- function(input, output, session) {
  
#Index File Display  
  
  #Inputting data
  reactindex <- reactive ({
    req(input$index)
    indexFile <- Biostrings::readDNAStringSet(input$index$datapath, format = "fastq")
    indexFile
  })
  #Updates slider
  observe({ 
    updateSliderInput (session,"displayindex",min=1,value=10,max=length(reactindex()))
  })
  #Display Index Data
  updatedisplayindex<- reactive({
    display.ind<-reactindex()[1:input$displayindex,]
    return(display.ind)
  })
  #Output Data
  output$tableindex <- renderTable({
    updatedisplayindex()[1:input$displayindex,]
  })

#Read File Display
  
  #Inputting data
  reactread <- reactive ({
    req(input$read)
    readFile<- Biostrings::readQualityScaledDNAStringSet(input$read$datapath)
    readFile
  })
  #Updates slider
  observe({ 
    updateSliderInput (session,"displayread",min=1,value=10,max=length(reactread()))
  })
  #Display Index Data
  updatedisplayread<- reactive({
    displayr<-reactread()[1:input$displayread,]
    return(displayr)
  })
  #Output Data
  output$tableread <- renderTable({
    updatedisplayread()[1:input$displayread,]
  })  
  
#Barcode File Display  
  
  #Inputting data
  reactbarcode <- reactive ({
    req(input$barcode)
    barcodeFile <- read.table(input$barcode$datapath, header = input$headerbarcode, sep = "\t")
    barcodeFile 
  })
  #Updates slider
  observe({ 
    updateSliderInput (session,"displaybarcode",min=1,value=10,max=nrow(reactbarcode()))
  })
  #Display Index Data
  updatedisplaybarcode<- reactive({
    displayb<-reactbarcode()[1:input$displaybarcode,]
    return(displayb)
  })
  #Output Data
  output$tablebarcode <- renderDataTable({
    updatedisplaybarcode()[1:input$displaybarcode,]
  })

#rcBarcodes
  
  rc<-reactive({
    return(as.logical(input$rcBarcodes))
 })  

#hammingDist
  
  ham<-reactive({
    return(as.numeric(input$hammingDist))
  })

#cores
  
 cores<-reactive({
    return(as.numeric(input$cores))
 })  
 
 out<-reactive({
    return(as.character(input$location))
 })  
 
#Demultiplex
observeEvent(input$demultiplex, {
  
  demultiplex(bcFile = reactbarcode() ,inds = reactindex() ,reads = reactread(),location=file.path(out(),"/Demultiplexed Sample Files"), rcBarcodes = rc(), cores= cores(), hammingDist= ham()) 
  
  })   
}  
shinyApp(ui, server)

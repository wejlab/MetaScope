
tabPanel(title = "Demultiplexing Read",
  sidebarLayout(
    sidebarPanel(
      #Read in Barcode File 
      
      fileInput("barcode", "Barcode File",
                accept = c(
                  ".txt")),
      checkboxInput("headerbarcode", "Header", TRUE),
      sliderInput("displaybarcode","View Rows:",min=0,max=0,value=10,step=1),
      
      #Read in Read and Index File          
      
      fileInput("index", "Index File",
                accept = c(
                  ".fast1.gz")),        
      fileInput("read", "Read File",
                accept = c(
                  ".fast1.gz")),
      sliderInput("display","View Rows:",min=0,max=0,value=10,step=1),
      actionButton("upload","Upload"),
      
      #Hamming, rcBarcode, and Codes
      
      selectInput("hammingDist","Hamming Distance",choices=1,selected=1),
      checkboxInput("rcBarcodes", "Reverse Compliment the Barcodes", FALSE),
      
      #Specify Output File
      
      textInput("location", "Output File Location",
                placeholder = "New Folder"
      ),
      
      #Dimultiplex!
      
      actionButton("demultiplex","Demultiplex")
    ),
    
    #Index Table Display 
    
    mainPanel(
      dataTableOutput("tablebarcode"), tableOutput("table")
    )
  )
)
)

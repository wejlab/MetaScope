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
           radioButtons("kingdom", "Choose a kingdom:",
                        c("Bacteria" = "bacteria",
                          "Virus" = "viral")
           ),
           
           # create checkbox input for representative library and reference library
           checkboxInput("representative", "representative", value = TRUE, width = NULL),
           checkboxInput("reference", "reference", value = FALSE, width = NULL),
           
           
           actionButton("downloadref","Download Ref_Seq")
         )
)

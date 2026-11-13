suppressPackageStartupMessages({
    library(shiny)
    library(ga.data)
  #  library(bettertrace)
    library(shinyjs)
    library(DT)
})

shinyUI(
  htmlTemplate(
    "www/data_processing.html",
    main_shiny_ui =
  fluidPage(

    useShinyjs(),
    ## Application title
    ## headerPanel("Plate Analyzer"),
    tags$div(
      style = "padding: 10px 15px 20px;",
      tags$button(
        id = "back_to_upload",
        type = "button",
        class = "circle-btn action-button",
        title = "Back to File Upload",
        HTML('<svg viewBox="0 0 24 24"><path d="M9 14L4 9l5-5"/><path d="M4 9h10.5a5.5 5.5 0 0 1 5.5 5.5 5.5 5.5 0 0 1-5.5 5.5H11"/></svg>')
      )
    ),
    tags$head(
        tags$style("#diPlot {height:60vh !important;}"),
        tags$link( rel="stylesheet", type="text/css", href="extra-styles.css" ),
        tags$style(HTML("
        .nav-group-container {
        position: fixed !important;
        top: 20px !important;
        
        /* Forces the container to ignore layout wrappers and stick to the true screen right */
        right: 20px !important;
        left: auto !important; 
        
        z-index: 99999 !important; /* Forces it over top header menus */
        
        display: flex !important;
        flex-direction: row !important; /* Keeps them horizontal */
        gap: 15px !important;           
        align-items: center !important;
        width: auto !important;         
        height: auto !important;
      }
      
      /* --- FLOATING CIRCLE CONTROL MECHANICS --- */
      .home-btn, .menu-btn {
        display: inline-flex !important;
        align-items: center !important;
        justify-content: center !important;
        width: 56px !important;
        height: 56px !important;
        border-radius: 50% !important;
        background-color: #ffffff !important;
        border: none !important;
        box-shadow: 5px 6px 14px rgba(0, 0, 0, 0.16) !important;
        cursor: default !important; 
        transition: all 0.2s ease-in-out !important;
        text-decoration: none !important;
        float: none !important;
        position: static !important;
        
        /* Reset any inherited float/position parameters from parent themes */
        float: none !important;
        position: static !important;
      }

      .home-btn:hover, .menu-btn:hover {
        transform: scale(1.04) !important;
        box-shadow: 6px 8px 18px rgba(0, 0, 0, 0.22) !important;
      }

      .home-icon-svg, .menu-icon-svg {
        width: 34px !important;
        height: 34px !important;
        fill: none !important;
        stroke: #22846C !important;
        stroke-width: 1.5 !important;
        stroke-linecap: round !important;
        stroke-linejoin: round !important;
      }
      
      .dropdown-wrapper {
        position: relative !important;
        display: inline-block !important;
      }

      .custom-dropdown-card {
        background: white !important;
        border-radius: 12px !important;
        box-shadow: 5px 6px 20px rgba(0, 0, 0, 0.16) !important;
        padding: 20px !important;
        width: 260px !important;
        position: absolute !important;
        
        /* Lock dropdown right edge perfectly under the hamburger button right edge */
        right: 0 !important;
        left: auto !important;
        top: 68px !important;
        
        text-align: left !important;
        z-index: 100000 !important;
      }
      
      .menu-row {
        display: flex !important;
        align-items: center !important;
        color: #22846C !important;
        font-size: 18px !important;
        font-weight: 500 !important;
        padding: 12px 8px !important;
        text-decoration: none !important;
        cursor: default !important; 
        border-radius: 6px !important;
      }
      .menu-row:hover { background-color: #f0fdfa !important; }
      .menu-row i { font-size: 22px !important; margin-right: 15px !important; }
      
      .menu-footer {
        display: flex !important;
        justify-content: center !important;
        gap: 25px !important;
        margin-top: 15px !important;
        padding-top: 15px !important;
        border-top: 1px solid #f3f4f6 !important;
        color: #9ca3af !important;
      }
      .footer-icon:hover { color: #4b5563 !important; cursor: pointer; }
      
      #uploaded_files_table + .shiny-input-container .progress {
        display: none !important;
      }
      
      .preview-container {
        width: 100%;
        max-width: 100%;
        box-sizing: border-box;
        overflow-x: auto;
        overflow-y: visible;
        padding: 0 5px;
        margin: 10px 0;
      }

/* Remove Bootstrap negative margins */
.preview-container .row {
    margin-left: 0 !important;
    margin-right: 0 !important;
}

/* Keep preview controls visible */
.preview-container .row > [class*=\"col-\"] {
    padding-left: 5px !important;
    padding-right: 5px !important;
}

/* Keep DT output within the panel */
.preview-container .dataTables_wrapper {
    width: 100% !important;
    max-width: 100% !important;
    box-sizing: border-box;
}
    .preview-toggle.preview-active,
    .preview-toggle.preview-active span {
     color: #198754 !important;
     font-weight: bold !important;
      }
    .file-name-status {
    display: inline-block;
    margin-left: 10px;
    color: #198754;
    font-weight: bold;
      }
      .modal-dialog.qc-wide-modal {
      width: 95% !important;
      max-width: 1600px !important;
      min-width: 90% !important;
      margin: 20px auto !important;
      }
     .qc-overview-btn {
      background-color: white !important;
      color: #444444 !important;
      font-size: 20px !important;
      font-weight: 700 !important;
      border: 2px solid #1D9685 !important;
      border-radius: 8px !important;
      padding: 16px 32px !important;
      min-width: 220px !important;
      min-height: 60px !important;
      cursor: default !important; 

      /* Shadow to the right and bottom */
      box-shadow: 5px 5px 8px rgba(0, 0, 0, 0.25) !important;
      transition: all 0.2s ease-in-out;
    }

    .qc-overview-btn:hover {
      background-color: #1D9685 !important;
      color: #f0f8f0 !important;
      border-color: #1D9685 !important;
      box-shadow: 7px 7px 12px rgba(0, 0, 0, 0.30) !important;
      transform: translate(-1px, -1px);
    }
    .qc-overview-btn:active {
      background-color: #444444 !important;
      box-shadow: 2px 2px 4px rgba(0, 0, 0, 0.20) !important;
      transform: translate(2px, 2px);
    }

    .qc-overview-btn:focus {
      color: #333333 !important;
      outline: none;
    }

    .modal-dialog.qc-wide-modal .modal-content {
      width: 100% !important;
      box-sizing: border-box !important;
    }

    .modal-dialog.qc-wide-modal .modal-body {
      overflow-x: auto !important;
      overflow-y: auto !important;
      max-height: 80vh !important;
    }
         .green-panel {
         background-color: #E8F5E9;
         border: 1px solid #A5D6A7;
         border-radius: 10px;
         padding: 20px;
         box-shadow: 4px 3px 8px rgba(0, 0, 0, 0.12);
         }
        .download-icon1 {
        display: inline-flex !important;
        flex-direction: row !important;
        align-items: center !important;
        justify-content: center !important;
        flex-wrap: nowrap !important;
        gap: 6px !important;

        background: transparent !important;
        border: none !important;
        padding: 2px 5px !important;
        color: #1D9685 !important;
        font-size: 14px !important;
        white-space: nowrap !important;
        width: auto !important;
        }
       .preview-toggle {
       background: transparent !important;
       border: none !important;
       box-shadow: none !important;
       color: #1D9685 !important;
       padding: 4px 8px !important;
       font-size: 16px !important;
       }
      .download-icon1 i {
      display: inline-block !important;
      margin: 0 !important;
      }

     .download-icon1 .download-label {
     display: inline-block !important;
     white-space: nowrap !important;
     }
     .download-icon2 {
     display: inline-flex !important;
     flex-direction: row !important;
     align-items: center !important;
     justify-content: center !important;
     flex-wrap: nowrap !important;
     gap: 6px !important;

     background: transparent !important;
     border: none !important;
     padding: 2px 5px !important;
     color: #FFA500 !important;
     font-size: 14px !important;
     white-space: nowrap !important;
     width: auto !important;
     }

     .download-icon2 i {
     display: inline-block !important;
     margin: 0 !important;
     }
    .download-icon2 .download-label {
    display: inline-block !important;
    white-space: nowrap !important;
    }
    .download-icon3 {
    display: inline-flex !important;
    flex-direction: row !important;
    align-items: center !important;
    justify-content: center !important;
    flex-wrap: nowrap !important;
    gap: 6px !important;
    background: transparent !important;
    border: none !important;
    padding: 2px 5px !important;
    color: #FFC0CB !important;
    font-size: 14px !important;
    white-space: nowrap !important;
    width: auto !important;
    }

   .download-icon3 i {
   display: inline-block !important;
   margin: 0 !important;
   }

  .download-icon3 .download-label {
  display: inline-block !important;
  white-space: nowrap !important;
  }
  .download-icon4 {
  display: inline-flex !important;
  flex-direction: row !important;
  align-items: center !important;
  justify-content: center !important;
  flex-wrap: nowrap !important;
  gap: 6px !important;
  background: transparent !important;
  border: none !important;
  padding: 2px 5px !important;
  color: #B8860B !important;
  font-size: 14px !important;
  white-space: nowrap !important;
  width: auto !important;
  }

 .download-icon4 i {
  display: inline-block !important;
  margin: 0 !important;
  }

 .download-icon4 .download-label {
  display: inline-block !important;
  white-space: nowrap !important;
 } 
 Shiny.addCustomMessageHandler('openPCPlot', function(message) {
    var win = window.open('', 'PCPlotWindow',
                          'width=1300,height=900,resizable=yes,scrollbars=yes');

    if (!win) {
      alert('Please allow pop-ups for this application.');
      return;
    }

    win.document.open();
    win.document.write(
      '<html><head><title>PC Plot</title></head>' +
      '<body style=\"font-family:Arial;padding:20px;\">' +
      '<h2>PC Plot</h2>' +
      '<p>Loading plot...</p>' +
      '<img id=\"pcplot-image\" style=\"max-width:100%;height:auto;\" />' +
      '</body></html>'
    );
    win.document.close();

    var img = win.document.getElementById('pcplot-image');
    img.src = message.src;

    img.onload = function() {
      var p = win.document.querySelector('p');
      if (p) p.remove();
    };
  });
 ")),
        tags$script(HTML("
    $(document).on('change', '#uploaded_files_table', function() {
      var files = this.files;
      var container = $(this).closest('.shiny-input-container');

      container.find('.file-name-status').remove();

      if (files.length > 0) {
        var names = [];

        for (var i = 0; i < files.length; i++) {
          names.push(
            '<div class=\"file-name-status\">&#10003; ' +
            $('<div>').text(files[i].name).html() +
            '</div>'
          );
        }

        $(this).closest('.form-group').append(names.join(''));
      }
    });
  "))
  ),
  
  # MENU BUTTON DROPDOWN WRAPPER
  div(class = "nav-group-container",
      
      # 1. Home Link
      tags$a(
        href = "https://gamap.genetic-analysis.com", 
        class = "home-btn",
        tags$svg(
          class = "home-icon-svg",
          viewBox = "0 0 24 24",
          tags$path(d = "M3 10.5 L12 3.5 L21 10.5 V19.5 C21 20.8 19.8 21 18.5 21 H5.5 C4.2 21 3 20.8 3 19.5 Z"),
          tags$path(d = "M10 21 V15.5 C10 14.2 11 13.5 12 13.5 C13 13.5 14 14.2 14 15.5 V21")
        )
      ),
      
      # 2. Hamburger Dropdown Action
      div(class = "dropdown-wrapper",
          tags$button(
            id = "btn_toggle",
            class = "menu-btn action-button", 
            tags$svg(
              class = "menu-icon-svg",
              viewBox = "0 0 24 24",
              tags$path(d = "M3 6h18"),
              tags$path(d = "M3 12h18"),
              tags$path(d = "M3 18h18")
            )
          ),
          
          # The Dropdown popout target
          conditionalPanel(
            condition = "input.btn_toggle % 2 == 1", 
            div(class = "custom-dropdown-card",
                
                # WIRED LINK 1: Marker Catalog (Quiet Excel download)
                tags$a(
                  href = "https://gamap.genetic-analysis.com/catalog.html",
                  class = "menu-row",
                  icon("table-cells-large"), 
                  "Marker catalog"
                ),
                
                # WIRED LINK 2: Documents 
                tags$a(
                  href = "https://gamap.genetic-analysis.com/discovery_documents.html", 
                  target = "_blank", 
                  class = "menu-row",
                  icon("file-shield"), 
                  "Documents"
                ),
                
                # Bottom decoration links
                div(class = "menu-footer",
                    icon("reply", class = "footer-icon"),
                    icon("sun", class = "footer-icon"),
                    icon("shield-check", class = "footer-icon")
                )
            )
          )
      )
  ),
  tags$div(
      style = paste(
        "display: flex;",
        "align-items: left;",
        "gap: 20px;",
        "padding: 5px 5px;"
      ),
  tags$div(
        style = "display:none;",
        plotOutput("PCplot_window", width = "1200px", height = "800px")
      ),
      
      tags$img(
        src = "gamap.png",
        style = "width: 480px; height: auto;"
      )
    ),

            tags$div(
              style = "display: none;",
              checkboxInput(
                "show_qc_page",
                label = NULL,
                value = FALSE
              )
            ),
            tags$div(
              style = "display: none;",
              checkboxInput(
                "show_pc_page",
                label = NULL,
                value = FALSE
              )
            ),
            conditionalPanel(
              condition = "!input.show_qc_page && !input.show_pc_page",
              
             wellPanel(
              div(
                class = "green-panel",
                fluidRow(

                    column(
                        3,
                        #fileInput("bc_file", label=h4("Select file/files"), multiple=TRUE,
                         #         accept="text/csv" )
                        h4("Uploaded Files"),
                        uiOutput("uploaded_files_list")
                    ),
                    column(4),
                    column(
                      2,
                      br(),
                      actionButton(
                        "qc_overview",
                        "QC Overview",
                        class = "qc-overview-btn"
                      )
                    ),
                    column(
                      3,
                      br(),
                      actionButton(
                        "pc_plot",
                        "PC Plot",
                        icon = icon("chart-scatter"),
                        class = "qc-overview-btn"
                      )
                    )
                    )
              )
            ),
            
            div(
              class = "preview-container",
              
              # ROW 1: Preview buttons and ALL download buttons
              fluidRow(
                shinyjs::hidden(
                column(
                  1,
                  actionButton(
                    "toggle_excel_preview1",
                    label = "Product 1",
                    icon = icon("eye-slash"),
                    class = "preview-toggle"
                  ),
                  id = "preview-button-column"
                )
                ),
                shinyjs::hidden(
                column(
                  1,
                  actionButton(
                    "toggle_excel_preview2",
                    label = "Product 2",
                    icon = icon("eye-slash"),
                    class = "preview-toggle"
                  ),
                  id = "preview-button-column2"
                )
                ),
                shinyjs::hidden(
                column(
                  1,
                  actionButton(
                    "toggle_excel_preview3",
                    label = "Product 3",
                    icon = icon("eye-slash"),
                    class = "preview-toggle"
                  ),
                  id = "preview-button-column3"
                )
                ),
                shinyjs::hidden(
                column(
                  1,
                  actionButton(
                    "toggle_excel_preview",
                    label = "Copy Number",
                    icon = icon("eye-slash"),
                    class = "preview-toggle"
                  ),
                  id = "preview-button-column4"
                )
                ),
                
                # Download all files
                shinyjs::hidden(
                  column(
                    1,
                    downloadButton(
                      "downloadResults",
                      label = "Download all files",
                      icon = icon("download"),
                      class = "download-icon1"
                    ),
                    id = "download-button-column"
                  )
                ),
                
                # Download Product 1
                shinyjs::hidden(
                  column(
                    1,
                    downloadButton(
                      "downloadResults2",
                      label = "Download Product 1",
                      icon = icon("download"),
                      class = "download-icon2"
                    ),
                    id = "download-button-column2"
                  )
                ),
                
                # Download Product 2
                shinyjs::hidden(
                  column(
                    1,
                    downloadButton(
                      "downloadResults3",
                      label = "Download Product 2",
                      icon = icon("download"),
                      class = "download-icon3"
                    ),
                    id = "download-button-column3"
                  )
                ),
                
                # Download Product 3
                shinyjs::hidden(
                  column(
                    1,
                    downloadButton(
                      "downloadResults4",
                      label = "Download Product 3",
                      icon = icon("download"),
                      class = "download-icon4"
                    ),
                    id = "download-button-column4"
                  )
                ),
                
                # Download Copy Number
                shinyjs::hidden(
                  column(
                    1,
                    downloadButton(
                      "downloadResults5",
                      label = "Download Copy Number",
                      icon = icon("download"),
                      class = "download-icon1"
                    ),
                    id = "download-button-column5"
                  )
                )
              ),
              
              # ROW 2: Preview tables (independent of download buttons)
              
              shinyjs::hidden(
                div(
                  id = "excel_preview_section1",
                  DT::DTOutput("product1")
                )
              ),
              
              shinyjs::hidden(
                div(
                  id = "excel_preview_section2",
                  DT::DTOutput("product2")
                )
              ),
              
              shinyjs::hidden(
                div(
                  id = "excel_preview_section3",
                  DT::DTOutput("product3")
                )
              ),
              
              shinyjs::hidden(
                div(
                  id = "excel_preview_section",
                  DT::DTOutput("excel_preview")
                )
              )
            )

        ),
        conditionalPanel(
          condition = "input.show_qc_page",
          
          fluidRow(
            column(4),
            column(
              4,
              actionButton(
                "qc_back",
                label = "Back to processed data",
                icon = icon("arrow-left"),
                class = "qc-overview-btn"
              ),
              actionButton(
                "qc_to_pc",
                "PC Plot",
                icon = icon("chart-scatter"),
                class = "qc-overview-btn"
              )
            )
          ),
          fluidRow(
            column(
              10,
              uiOutput("qcProductsAll")
            )
          ),
          
          fluidRow(
            column(
              10,
              DT::dataTableOutput("qcTableST")
            )
          ),
          
          fluidRow(
            column(
              10,
              DT::dataTableOutput("qcTableQCC29")
            )
          ),
          
          fluidRow(
            column(
              10,
              DT::dataTableOutput("qcTableQCC33")
            )
          ),
          fluidRow(
            column(
              10,
              DT::dataTableOutput("qcTableZymo")
            )
          )
        ),
        conditionalPanel(
          condition = "input.show_pc_page",
          
          fluidRow(
            column(4),
            column(
              4,
              actionButton(
                "pc_back",
                label = "Back to processed data",
                icon = icon("arrow-left"),
                class = "qc-overview-btn"
              ),
              actionButton(
                "pc_to_qc",
                "QC Overview",
                icon = icon("arrow-left"),
                class = "qc-overview-btn"
              )
            )
          ),
          
          fluidRow(
            column(
              10,
              plotOutput("PCplot", height = "700px")
            )
          )
        )
       
 )
  )
)

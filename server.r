options( stringsAsFactors=FALSE )

suppressPackageStartupMessages({
    library(shiny)
    library(shinyjs)
    library(devtools)
    library(ga.data)
    library(ga.gamap)
    library(ga.utils)
    library(ga.gamapqc)
    library(foreach)
    library(ggplot2)
    library(purrr)
    library(stringr)
    library(ga.software)
    library(ga.software.dd)
    library(ga.biocode)
    library(ga.luminex)
    library(dplyr)
    library(tibble)
    library(ga.biobank)
    library(cowplot)
    library(ga.discovery)
    #library(tidyverse)
})



colMedians <- function(x, ...){
  apply(x, 2, median, ...)
}

probe.ibd.re <- function (include.technical = FALSE, anchor = FALSE) {
  rx <- "(DP\\d{4})"
  if (include.technical) 
    rx <- "(DP\\d{4}|HYC\\d+|HYPCON+|UL+)"
  if (anchor) 
    rx <- paste0("^", rx, "$")
  rx
}

ibd.bacteria <- function(probelist, ibd.probe){
  if(  any( probelist$Probe %in% ibd.probe ) ){
    bacteria <- probelist[ibd.probe, "Bacteria"]
    return( bacteria)
  } else { print("The probe is not in the list") }
}

#PC1 and PC2:
plot.pca.scores.12 <- function(pca.obj,d,title) {
  ggplot( cbind(d,vegan::scores(pca.obj)[]), aes(x=PC1,y=PC2) ) +
    pca.axis.labels(pca.obj) +
    geom_point( aes( col = Sample) ) +
    theme_classic() +
    ggtitle( paste("PC1,PC2; ", title)) +
    geom_blank()
}


DEBUG <- FALSE

## Define server logic required to draw a histogram
shinyServer(function(input, output, session) {



    input_file <- function() {
        b <- input$bc_file
        
        b
    }
    # work on function that will merge all probes
    input_folder <- function() {
        
        b <- dir(input$bc_file, full.names = TRUE)
        
    }


    ## FUNCTIONS

   plateData<-reactive({

        bc.file <- req(input_file())

        x_run=discovery(bc.file$datapath,stop.at="file")
       x_run
    })
    # function that will merge all probes
    plateDataFromFolder<-reactive({

        bc.file <- req(input_folder())

        x_run=read.luminex(bc.file$datapath)
       x_run
    })

    qctable <- function(qn) {
      
      pd <- req(plateData())
      
      pd2 <- pd[ grep( paste0("^",qn), pd$Sample ),  ]
      
      i <- grep( "HYPCON|HYC01|UL", colnames(pd2), value=TRUE )
      pd3 <- pd2[, c("Sample",i)]
    
      
      
      
      #pd3$Total <- sprintf( "%.0f", rowSums(xr, na.rm=TRUE) )
      
      o <- order(
        pd3$Sample %!~% "^ST",
        pd3$Sample %!~% "^QCC29",
        pd3$Sample %!~% "^QCC23",
        pd3$Sample %!~% "^QCC33",
        pd3$Sample %!~% "^Zymo"
      )
      
      pd4 <- pd3[o,]
      
     
      
      cn <- colnames(pd4)
      
      pd5 <- foreach( i=1:ncol(pd4), .combine = cbind.data.frame ) %do% {
        v <- pd4[,i]
        d <- data.frame(v=v)
        colnames(d) <- cn[i]
        
        if( cn[i] == "HYPCON|HYC01" ) {
          d <- cbind.data.frame( d, "" )
          colnames(d)[2] <- "?"
        }
        
        
        d
        
      }
      
      colnames(pd5)[ colnames(pd5) == "?" ] <- ""
      pd5
      
    }
    
    normalized<-reactive({
      Set1x1 <- req(plateData())
      Set1x <- as.matrix( Set1x1[ , grepl( probe.ibd.re(include.technical=TRUE), colnames(Set1x1) ) ] )      
      rownames(Set1x) <- Set1x1$Sample
     
      
      ## HYC-normalized data:
      Set1.hyc.i <-discovery(Set1x1,start.from="file",stop.at="hyb",qc.check.UL=FALSE)
       
      ## HYC and ST5-normalized data:
      Set1.st5.i<-discovery(Set1x1,start.from="file",stop.at="qcc",qc.check.UL=FALSE)
      
      # correct for QCC29 mean
      Set1.i.st5.qcc29<-discovery(Set1x1,start.from="file",stop.at="qcc29",qc.check.UL=FALSE)
      rownames(Set1.i.st5.qcc29) <- Set1x1$Sample
      Set1.i.st5.qcc29
      
    })
    # PLOT EACH PROBE RAW SIGNAL FOR EACH GROUP
    output$Discoveryplot<-renderPlot({ 

    Set1x1<-req(plateData())
    Set1.raw<-NULL
   
    Set1x <- as.matrix( Set1x1[ , grepl( probe.ibd.re(include.technical=TRUE), colnames(Set1x1) ) ] )
    rownames(Set1x) <- Set1x1$Sample
    idx<-Set1x[,"UL"]<50
    Set1.raw <- rbind(Set1.raw,Set1x)
    Set1.df <- cbind.data.frame(Set1x1$Sample,Set1.raw[, order(colnames(Set1.raw), decreasing = FALSE)])
    Set1.dfm <- reshape2::melt(Set1.df)
    colnames(Set1.dfm) <- c("Sample","Probe", "Raw_signal")
    ist<-grepl("ST|^Zymo_FR|QCC23|QCC33",Set1.dfm$Sample)
    Set1.dfm<-Set1.dfm[ist,]
    # put all multiple controls to one:
    Set1.dfm$Sample<-gsub("_A","",Set1.dfm$Sample)
    Set1.dfm$Sample<-gsub("_B","",Set1.dfm$Sample)
    Set1.dfm$Sample<-gsub("_a","",Set1.dfm$Sample)
    Set1.dfm$Sample<-gsub("_b","",Set1.dfm$Sample)
    Set1.dfm$Sample<-gsub("_c","",Set1.dfm$Sample)
    g1=ggplot(Set1.dfm, aes( x = Probe, y = Raw_signal) ) +
      geom_point(aes(col=Sample))+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7))+
    # geom_jitter( aes( col = Group, shape = Group), position = position_jitter(w = 0.2, h = 0.2)) +
    # scale_shape_manual(values = 1:length(unique(Set1.dfm$Group))) +
     ggtitle("Probe set 1 - Raw signals") 
    # coord_flip()
    #g2=ggplot(Set1.dfm,aes( x=Probe, y=Raw_signal,color=Group))  + geom_point()
    #plot_grid(g1,g2,labels="AUTO",align="h",ncol=1)
    print(g1)
    })
     
    output$downloadResults <- downloadHandler(
      
      filename = function() {
        sub( "\\.csv$", "-normalized.csv", input_file()$name, ignore.case=TRUE )
      },
      
      content = function(file) {
        d <- req(normalized())
        d <- d[ !grepl("^ST|^QCC29|^QCC23|^QCC33|^Zymo",rownames(d)),] 
        dd <- cbind.data.frame( Sample=rownames(d), d)
        write.csv( dd, file, row.names=FALSE )
      }
      
    )
    
    observeEvent( input_file(), {
      if( is.null(input_file()) ) {
        shinyjs::hideElement( id="download-button-column" )
      } else {
        shinyjs::showElement( id="download-button-column" )
      }
    })  
    
   # PAGE 2
    ## Page 2: QC tables
    
    output$qcTableST <- renderTable(qctable("^ST"), digits=1)
    output$qcTableQCC29 <- renderTable(qctable("^QCC29"), digits=1)
    output$qcTableQCC23 <- renderTable(qctable("^QCC23"), digits=1)
    output$qcTableQCC33 <- renderTable(qctable("^QCC33"), digits=1)
    output$qcTableZymo <- renderTable(qctable("^Zymo_FR"), digits=1)
    
   # PAGE 3: PLOT PC PLOT FOR ALL SAMPLES: NORMALIZED DATA
    output$PCplot<-renderPlot({ 
      
      Set1x1<-req(plateData())
      Set1x <- as.matrix( Set1x1[ , grepl( probe.ibd.re(include.technical=TRUE), colnames(Set1x1) ) ] )
      rownames(Set1x) <- Set1x1$Sample
      idx<-Set1x[,"UL"]<50
      i.qcc<-grepl("QCC23|QCC33|ST|Zymo|QCC29",Set1x1$Sample)
      data<-req(normalized())
      data<-data[idx | !i.qcc,]
      data[data < 5] <- 2.5
     
      signal.data <- as.data.frame(data)
      pr.sd<-apply(signal.data,2,sd)
      i.sd0<-pr.sd==0
      signal.data<-signal.data[,!i.sd0]
      d<-as.data.frame(Set1x1$Sample[idx | !i.qcc])
      colnames(d)<-"Sample"
      pca <- prcomp(log2(signal.data))
      t <- "log data"
      
      p=plot.pca.scores.12(pca,d,title = t)
      p
    })
    
    
    
    ## Page 4: Probe Annotations

    load("C:/b/Projects/Development/Discovery/Data/ProcessedData/Probelist_new.Rmd")

    probe.a <- data.frame(
        Probe = probelist_new$Probe.name,
        Annotation = probelist_new$Probe.annotation,
        Genus = probelist_new$Genus,
        Family = probelist_new$Family,
        Order = probelist_new$Order,
        Class = probelist_new$Class,
        Phylum = probelist_new$Phylum
        
    )

    output$ProbeAnnotations <- renderDataTable(
        probe.a,
        options = list(
            pageLength = 50
        )
    )

})


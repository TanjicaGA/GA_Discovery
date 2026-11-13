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
    library(writexl)
    library(DT)
    library(tidyverse)
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

probe.ibd.re2 <- function (include.technical = FALSE, anchor = FALSE) {
  rx <- "(DP\\d{4})"
  if (include.technical) 
    rx <- "(DP\\d{4}|HYC\\d+|HYPCON+)"
  if (anchor) 
    rx <- paste0("^", rx, "$")
  rx
}
probe.ibd.re3 <- function (include.technical = FALSE, anchor = FALSE) {
  rx <- "(DP\\d{4})"
  if (include.technical) 
    rx <- "(DP\\d{4}|HYC\\d+)"
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
    ggtitle( paste("PC1,PC2 for all files; ", title)) +
    geom_blank()
}


DEBUG <- FALSE

## Define server logic required to draw a histogram
shinyServer(function(input, output, session) {



    input_file <- function() {
        b <- input$bc_file
        print("Someone is uploading files!!!!")        
        b
    }
    

    ## FUNCTIONS
    
    plateData0<-reactive({
      if(length(input$bc_file[,1])>=1){
      
      for(nr in 1:length(input$bc_file[,1])){
        
        data<-input$bc_file[[nr,'datapath']]
        Set1x1<-discovery(data,stop.at="file")
        Set1x1$Sample<-gsub("^ST\\d|^ST\\d\\d|^SynCon\\d\\d|SynCon\\d|ST$","ST",as.character(Set1x1$Sample))
        Set1x1$Sample<-gsub("^QCC29\\d|^QCC29\\d\\d|QCC29$","QCC29",as.character(Set1x1$Sample))
        Set1x1$Sample<-gsub("^QCC33\\d|^QCC33\\d\\d|QCC33$","QCC33",as.character(Set1x1$Sample))
        Set1x1$Sample<-gsub("^Zymo\\d|^Zymo\\d\\d|Zymo$","Zymo",as.character(Set1x1$Sample))
        if ("DP0049" %in% colnames(Set1x1)){
          colnames(Set1x1)[which(names(Set1x1)=="HYPCON")]<-"HYPCON_1"
          colnames(Set1x1)[which(names(Set1x1)=="UL")]<-"UL_1"
          colnames(Set1x1)[which(names(Set1x1)=="HYC01")]<-"HYC01_1"
        }
        if ("DP0004" %in% colnames(Set1x1)){
          colnames(Set1x1)[which(names(Set1x1)=="HYPCON")]<-"HYPCON_2"
          colnames(Set1x1)[which(names(Set1x1)=="UL")]<-"UL_2"
          colnames(Set1x1)[which(names(Set1x1)=="HYC01")]<-"HYC01_2"
        }
        if ("DP0200" %in% colnames(Set1x1)){
          
          colnames(Set1x1)[which(names(Set1x1)=="HYPCON")]<-"HYPCON_3"
          colnames(Set1x1)[which(names(Set1x1)=="UL")]<-"UL_3"
          colnames(Set1x1)[which(names(Set1x1)=="HYC01")]<-"HYC01_3"
          
        }
       
        if(nr==1){
          signal.data<-Set1x1
        }
        else if(nr>1){
          signal.data<-merge(signal.data,Set1x1,by="Sample")
        }
      }
      
      signal.data
      }
      
    })  

   hycData<-reactive({
    if(length(input$bc_file[,1])>=1){
     for(nr in 1:length(input$bc_file[,1])){
       
       data<-input$bc_file[[nr,'datapath']]
       Set1x1<-discovery(data,stop.at="file")
       Set1x1$Sample<-gsub("^ST\\d|^ST\\d\\d|^SynCon\\d\\d|SynCon\\d","ST",as.character(Set1x1$Sample))
       Set1.hyc.i <-discovery(Set1x1,start.from="file",stop.at="hyb",qc.check.UL=FALSE)
       rownames(Set1.hyc.i) <- Set1x1$Sample
       i.qcc<-grepl("QCC23|QCC33|ST|Zymo|QCC29",Set1x1$Sample)
       #Set1.hyc.i<-Set1.hyc.i[idx|!i.qcc,]
       Set1.hyc.i<-Set1.hyc.i[!i.qcc,]
       if(dim(Set1.hyc.i)[1]==0){stop(paste("No patient samples in the file!"))}
       Set1.hyc.df<-as.data.frame(Set1.hyc.i) %>% rownames_to_column("Sample")
       
       if(nr==1){
         signal.data<-Set1.hyc.df
       }
       else{
         signal.data<-merge(signal.data,Set1.hyc.df,by="Sample",sort=FALSE)
       }
     }
     sorted_columns<-c("Sample",sort(colnames(signal.data)[2:length(colnames(signal.data))]))
     signal.data<-signal.data[,sorted_columns]
    }
     
    })
   
   interpolationToCopyNumber<-reactive({
   d_df<-req(hycData())
   d_df$UL<-NULL
   d_df$HYC01<-NULL
   d_df$HYPCON<-NULL
   d_df$UL.x<-NULL
   d_df$HYC01.x<-NULL
   d_df$HYPCON.x<-NULL
   d_df$UL.y<-NULL
   d_df$HYC01.y<-NULL
   d_df$HYPCON.y<-NULL
   
   
   # add specificities:
   DF <- read.csv("C:/b/Projects/Development/Discovery/Data/ProcessedData/DIS06-04_Appendix_R1/probe_level_results.csv")
   interpolation_probes<-colnames(d_df)
   interpolation_probes<-interpolation_probes[startsWith(interpolation_probes,"DP")]
   bf <- bind_rows(lapply(sort(interpolation_probes), function(pr) {
     q <- DF |> filter(stage == "PCR", probe == pr, dilution %in% c(1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9)) |>
       arrange(copies)
     if(any(d_df[,pr]>min(q$median_signal)) & any(d_df[,pr])<max(q$median_signal)){
       tibble(probe = pr, 
              samples = d_df$Sample,
              estimated_copies_loglog = exp(approx(x = log(q$median_signal), y = log(q$copies), xout = log(d_df[,pr]))$y))
     }
     else{
       tibble(probe = pr, 
              samples = d_df$Sample,
              estimated_copies_loglog = rep(NA,length(rownames(d_df))))
     }
   }))
   bfreshaped <- bf %>% pivot_wider(
     names_from = `probe`,
     values_from = `estimated_copies_loglog`
   )
   })
   plateDataprobelist_set3<-reactive({
   d_df<-req(plateData())
   d_df$UL<-NULL
   d_df$HYC01<-NULL
   d_df$HYPCON<-NULL
   d_df$UL.x<-NULL
   d_df$HYC01.x<-NULL
   d_df$HYPCON.x<-NULL
   d_df$UL.y<-NULL
   d_df$HYC01.y<-NULL
   d_df$HYPCON.y<-NULL
   # add probe names as well:
   load("probelist_set3.Rmd")
   d_df <- d_df[,c("Sample",intersect(colnames(d_df),probelist_set3$Probe.name))]
   bn<-look.up( colnames(d_df)[2:length(d_df)], probelist_set3, to="Probe.annotation", leave.unfound=FALSE ) 
   bf<-rbind(c("Sample",bn),d_df)
   })
   plateDataprobelist_set2<-reactive({
   d_df<-req(plateData())
   d_df$UL<-NULL
   d_df$HYC01<-NULL
   d_df$HYPCON<-NULL
   d_df$UL.x<-NULL
   d_df$HYC01.x<-NULL
   d_df$HYPCON.x<-NULL
   d_df$UL.y<-NULL
   d_df$HYC01.y<-NULL
   d_df$HYPCON.y<-NULL
   # add probe names as well:
   load("probelist_set2.Rmd")
   d_df <- d_df[,c("Sample",intersect(colnames(d_df),probelist_set2$Probe.name))]
   bn<-look.up( colnames(d_df)[2:length(d_df)], probelist_set2, to="Probe.annotation", leave.unfound=FALSE ) 
   bf<-rbind(c("Sample",bn),d_df)
   })
   plateDataprobelist_set1<-reactive({
     d_df<-req(plateData())
     d_df$UL<-NULL
     d_df$HYC01<-NULL
     d_df$HYPCON<-NULL
     d_df$UL.x<-NULL
     d_df$HYC01.x<-NULL
     d_df$HYPCON.x<-NULL
     d_df$UL.y<-NULL
     d_df$HYC01.y<-NULL
     d_df$HYPCON.y<-NULL
     # add probe names as well:
     load("probelist_set1.Rmd")
     d_df <- d_df[,c("Sample",intersect(colnames(d_df),probelist_set1$Probe.name))]
     bn<-look.up( colnames(d_df)[2:length(d_df)], probelist_set1, to="Probe.annotation", leave.unfound=FALSE ) 
     bf<-rbind(c("Sample",bn),d_df)
   })
   plateData<-reactive({
     if(length(input$bc_file[,1])>=1){
       for(nr in 1:length(input$bc_file[,1])){
         
         data<-input$bc_file[[nr,'datapath']]
         Set1x1<-discovery(data,stop.at="file")
         Set1x1$Sample<-gsub("^ST\\d|^ST\\d\\d|^SynCon\\d\\d|SynCon\\d","ST",as.character(Set1x1$Sample))
         Set1.hyc.i <-discovery(Set1x1,start.from="file",stop.at="qcc29",qc.check.UL=FALSE)
         rownames(Set1.hyc.i) <- Set1x1$Sample
         i.qcc<-grepl("QCC23|QCC33|ST|Zymo|QCC29",Set1x1$Sample)
         #Set1.hyc.i<-Set1.hyc.i[idx|!i.qcc,]
         Set1.hyc.i<-Set1.hyc.i[!i.qcc,]
         if(dim(Set1.hyc.i)[1]==0){stop(paste("No patient samples in the file!"))}
         Set1.hyc.df<-as.data.frame(Set1.hyc.i) %>% rownames_to_column("Sample")
         
         if(nr==1){
           signal.data<-Set1.hyc.df
         }
         else{
           signal.data<-merge(signal.data,Set1.hyc.df,by="Sample",sort=FALSE)
         }
       }
       sorted_columns<-c("Sample",sort(colnames(signal.data)[2:length(colnames(signal.data))]))
       signal.data<-signal.data[,sorted_columns]
     }
     
   })
   


  qctable <- function(qn) {
      
      pd <- req(plateData0())
      
      pd2 <- pd[ grep( paste0("^",qn), pd$Sample ),  ]
      if(dim(pd2)[1]==0){
        pd2 <- pd[ grep(qn, pd$Sample ),  ]
      }
      if(dim(pd2)[1]!=0){
        # start with a logical statement:
        pass<-TRUE
      
        i <- grep( "HYPCON|HYC01|UL", colnames(pd2), value=TRUE )
        pd3 <- pd2[, c("Sample",i)]
        print(pd2)
        ### total QCC33, total QCC29, total ST:
        if(qn=="QCC33"){
          for(nr in 1:length(input$bc_file[,1])){
            data<-input$bc_file[[nr,'datapath']]
            Set1x1<-discovery(data,stop.at="file")
            Set_qcc33 <- as.matrix(Set1x1[,grepl( probe.ibd.re2(include.technical=TRUE), colnames(Set1x1) ) ] )
            rownames(Set_qcc33) <- Set1x1$Sample
            ## QCC33 data:
           
            pd_qcc33 <- Set_qcc33[ grep( paste0("^",qn), rownames(Set_qcc33) ),  ]
            if(length(pd_qcc33)==0){
              pd_qcc33 <- Set_qcc33[ grep(qn, rownames(Set_qcc33) ),  ]
            }
            if(is.null(dim(pd_qcc33))){
              qcc33_total<-sum(pd_qcc33)
            }
            else{
              qcc33_total<-rowSums(pd_qcc33, na.rm=TRUE)
            }
            if ("DP0049" %in% colnames(Set1x1)){ 
              pass_1<-TRUE
              for(i in 1:length(qcc33_total)){
                 if(qcc33_total[i]>40000 && qcc33_total[i]<68000){
                   pass<-pass & TRUE
                   pass_1<-pass_1 & TRUE
                 }
                 else{
                   pass<-pass & FALSE
                   pass_1<-pass_1 & FALSE
                
                 }
              }
              pd3$`QCC33 Total_1`<-qcc33_total
              
            }
            if ("DP0004" %in% colnames(Set1x1)){
              pass_2<-TRUE
              for(i in 1:length(qcc33_total)){
                if(qcc33_total[i]>15000 && qcc33_total[i]<52000){
                  pass<-pass & TRUE
                  pass_2<-pass_2 & TRUE
                }
                else{
                  pass<-pass & FALSE
                  pass_2<-pass_2 & FALSE
                
                }
              }
              pd3$`QCC33 Total_2`<-qcc33_total
            }
            if ("DP0200" %in% colnames(Set1x1)){
              pass_3<-TRUE
              for(i in 1:length(qcc33_total)){
                if(qcc33_total[i]>9000 && qcc33_total[i]<100100){
                  pass<-pass & TRUE
                  pass_3<-pass_3 & TRUE
                }
                else{
                  pass<-pass & FALSE
                  pass_3<-pass_3 & FALSE
                
                }
              }
              pd3$`QCC33 Total_3`<-qcc33_total
             
            }
          }
        }
        if(qn=="ST"){ 
          for(nr in 1:length(input$bc_file[,1])){
        
            data<-input$bc_file[[nr,'datapath']]
            Set1x1<-discovery(data,stop.at="file")
            Set_qcc30 <- as.matrix(Set1x1[,grepl( probe.ibd.re2(include.technical=TRUE), colnames(Set1x1) ) ] )
            rownames(Set_qcc30) <- Set1x1$Sample
            
            ## ST5 data:
            pd_qcc30 <- Set_qcc30[ grep( paste0("^",qn), rownames(Set_qcc30) ),  ]
            if(length(pd_qcc30)==0){
              pd_qcc30 <- Set_qcc30[ grep(qn, rownames(Set_qcc30) ),  ]
            }
            if(is.null(dim(pd_qcc30))){
               st_total<-sum(pd_qcc30)
            }
            else{
               st_total<-rowSums(pd_qcc30, na.rm=TRUE)
            }   
            if ("DP0049" %in% colnames(Set1x1)){
              pass_1<-TRUE
              for(i in 1:length(st_total)){
                if(st_total[i]>530000 && st_total[i]<660000){
                  pass<-pass & TRUE
                  pass_1<-pass_1 & TRUE
                 }
                else{
                  pass<-pass & FALSE
                  pass_1<-pass_1 & FALSE
            
                 }
              }
              pd3$`ST Total_1`<-st_total
             }
            if ("DP0004" %in% colnames(Set1x1)){
              pass_2<-TRUE
              for(i in 1:length(st_total)){
                if(st_total[i]>420000 && st_total[i]<550000){
                  pass<-pass & TRUE
                  pass_2<-pass_2 & TRUE
                }
                else{
                  pass<-pass & FALSE
                  pass_2<-pass_2 & FALSE
            
                }
              }
              pd3$`ST Total_2`<-st_total
             }
            if ("DP0200" %in% colnames(Set1x1)){
              pass_3<-TRUE
              for(i in 1:length(st_total)){
               if(st_total[i]>380000 && st_total[i]<500000){
                 pass<-pass & TRUE
                 pass_3<-pass_3 & TRUE
               }
               else{
                 pass<-pass & FALSE
                 pass_3<-pass_3 & FALSE
            
               }
              }
             pd3$`ST Total_3`<-st_total
            }
          }
        }
        if(qn=="QCC29"){
          for(nr in 1:length(input$bc_file[,1])){
          
            data<-input$bc_file[[nr,'datapath']]
            Set1x1<-discovery(data,stop.at="file")
            Set_qcc29 <- as.matrix(Set1x1[,grepl( probe.ibd.re3(include.technical=TRUE), colnames(Set1x1) ) ] ) 
            rownames(Set_qcc29) <- Set1x1$Sample
            ## QCC29 data:
            #i.qcc29 <- grepl("^QCC29", rownames(Set_qcc29), ignore.case = TRUE)
            #qcc29_total<-sum(colMedians( Set_qcc29[i.qcc29 , ,drop=FALSE ], na.rm=TRUE ))
            pd_qcc29 <- Set_qcc29[ grep( paste0("^",qn), rownames(Set_qcc29) ),  ]
            if(length(pd_qcc29)==0){
              pd_qcc29 <- Set_qcc29[ grep(qn, rownames(Set_qcc29) ),  ]
            }
            if(is.null(dim(pd_qcc29))){
               qcc29_total<-sum(pd_qcc29) 
            }
            else{
               qcc29_total<-rowSums(pd_qcc29, na.rm=TRUE)
            }
            if("DP0049" %in% colnames(Set1x1)){
              pd3$`QCC29 Total_1`<-qcc29_total
              pass_1<-TRUE
              if(all(qcc29_total<10000)){
                pass_1<-pass_1 & TRUE
              }
              else{
                pass_1<-pass_1 & FALSE
                
              }
            }
            if("DP0004" %in% colnames(Set1x1)){
              pd3$`QCC29 Total_2`<-qcc29_total
              pass_2<-TRUE
              if(all(qcc29_total<10000)){
                pass_2<-pass_2 & TRUE
              }
              else{
                pass_2<-pass_2 & FALSE
                
              }
            }
            if("DP0200" %in% colnames(Set1x1)){
              pd3$`QCC29 Total_3`<-qcc29_total
              pass_3<-TRUE
              if(all(qcc29_total<10000)){
                pass_3<-pass_3 & TRUE
              }
              else{
                pass_3<-pass_3 & FALSE
                
              }
            }
            if(all(qcc29_total<10000)){
              pass<-pass & TRUE
            }
            else{
              pass<-pass & FALSE
          
           }
        
         }
       }
      
        if (qn=="QCC29"){
         if ("HYC01_1" %in% colnames(pd) ){ 
           pass_1==TRUE
          if((all(pd2[,"HYC01_1"]>1500) && all(pd2[,"HYC01_1"]<3200) && all(pd2[,"UL_1"]<50) && all(pd2[,"HYPCON_1"]>400) && all(pd2[,"HYPCON_1"]<1250))){
            pass<-pass & TRUE
            pass_1<-pass_1 & TRUE
          }  
          else{
            pass<-pass & FALSE
            pass_1<-pass_1 & FALSE
         
          }
          }
         if ("HYC01_2" %in% colnames(pd) ){   
           pass_2==TRUE
          if((all(pd2[,"HYC01_2"]>1500) && all(pd2[,"HYC01_2"]<3200) && all(pd2[,"UL_2"]<50) && all(pd2[,"HYPCON_2"]>400) && all(pd2[,"HYPCON_2"]<1250))){
            pass<-pass & TRUE
            pass_2<-pass_2 & TRUE
          }  
          else{
            pass<-pass & FALSE
            pass_2<-pass_2 & FALSE
           
          }
         }
         if ("HYC01_3" %in% colnames(pd) ){
           pass_3==TRUE
          if((all(pd2[,"HYC01_3"]>1500) && all(pd2[,"HYC01_3"]<3200) && all(pd2[,"UL_3"]<50) && all(pd2[,"HYPCON_3"]>400) && all(pd2[,"HYPCON_3"]<1250))){
            pass<-pass & TRUE
            pass_3<-pass_3 & TRUE
          }
          else{
            pass<-pass & FALSE
            pass_3<-pass_3 & FALSE
           
          }
         }
        
        }
        else{
         if ("HYC01_1" %in% colnames(pd) ){  
           pass_1==TRUE
          if((all(pd2[,"HYC01_1"]>1500) && all(pd2[,"HYC01_1"]<3200) && all(pd2[,"UL_1"]<50) && all(pd2[,"HYPCON_1"]>4000))){
           pass<-pass & TRUE
           pass_1<-pass_1 & TRUE
          } 
          else{
           pass<-pass & FALSE
           pass_1<-pass_1 & FALSE
         }
         }
         if ("HYC01_2" %in% colnames(pd) ){ 
           pass_2==TRUE
          if((all(pd2[,"HYC01_2"]>1500) && all(pd2[,"HYC01_2"]<3200) && all(pd2[,"UL_2"]<50) && all(pd2[,"HYPCON_2"]>4000))){
           pass<-pass & TRUE
           pass_2<-pass_2 & TRUE
          }  
          else{
           pass<-pass & FALSE
           pass_2<-pass_2 & FALSE
          }
         }
         if ("HYC01_3" %in% colnames(pd) ){
           pass_3==TRUE
          if((all(pd2[,"HYC01_3"]>1500) && all(pd2[,"HYC01_3"]<3200) && all(pd2[,"UL_3"]<50) && all(pd2[,"HYPCON_3"]>4000))){
           pass<-pass & TRUE
           pass_3<-pass_3 & TRUE
           }
          else{
           pass<-pass & FALSE
           pass_3<-pass_3 & FALSE
           }
         }
        } 
        if(exists("pass_1")){
          if(pass_1==TRUE){
            pd3$Result_1<-"Pass"
          }
          else{
            pd3$Result_1<-"Fail"
          }
          
        }
        if(exists("pass_2")){
          if(pass_2==TRUE){
            pd3$Result_2<-"Pass"
          }
          else{
            pd3$Result_2<-"Fail"
          }
          
        }
        if(exists("pass_3")){
          if(pass_3==TRUE){
            pd3$Result_3<-"Pass"
          }
          else{
            pd3$Result_3<-"Fail"
          }
          
        }
      
        
       if(pass==TRUE){
         pd3$Result<-"Pass"
       }
       else{
        pd3$Result<-"Fail"
       }
       cn <- colnames(pd3)
       pd5 <- foreach( i=1:ncol(pd3), .combine = cbind.data.frame ) %do% {
       v <- pd3[,i]
       d <- data.frame(v=v)
       colnames(d) <- cn[i]
        
       if( cn[i] == "HYPCON|HYC01" ) {
          d <- cbind.data.frame( d, "" )
          colnames(d)[2] <- "?"
       }
        
        d
       }
      
      colnames(pd5)[ colnames(pd5) == "?" ] <- ""
      cols_1 <- cn[grepl("_1$", cn)]
      cols_2 <- cn[grepl("_2$", cn)]
      cols_3 <- cn[grepl("_3$", cn)]
      
      
      pd5 <- pd5[, c("Sample",cols_1, cols_2, cols_3,"Result"), drop = FALSE]
      
      
      
      
      # Preserve the column names as a separate row
      header2 <- names(pd5)
      
      # Create display data frame
      pd5_display <- rbind(
        header2,
        as.data.frame(lapply(pd5, as.character), check.names = FALSE)
      )
      pd5_display
    }
 }
    
    #normalized<-reactive({
    #  Set1x1 <- req(plateData())
      
      # correct for QCC29 mean
    #  Set1.i.st5.qcc29<-discovery(Set1x1,start.from="file",stop.at="qcc29",qc.check.UL=FALSE)
    #  rownames(Set1.i.st5.qcc29) <- Set1x1$Sample
    #  Set1.i.st5.qcc29
      
    #})
    # PLOT EACH PROBE RAW SIGNAL FOR EACH GROUP
    output$Discoveryplot<-renderPlot({ 

   
    Set1.df<-req(plateData())
    
    Set1.dfm <- reshape2::melt(Set1.df)
    colnames(Set1.dfm) <- c("Sample","Probe", "Normalized_signal")
    #remove low signals
    Set1.dfm<-Set1.dfm[Set1.dfm$Normalized_signal>100,]
    
    g1=ggplot(Set1.dfm, aes( x = Probe, y = Normalized_signal) ) +
     
      geom_point(aes(col=Sample),show.legend = FALSE)+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7))+
      
     ggtitle("Normalized signals higher than 100 for all plates") 
    
    print(g1)
    if(any(is.na(Set1.df))==TRUE){stop(paste("File is corrupted, some fields are empty.","\n",
                                       "However, you can still download normalized data."))}
    })
     
    output$downloadResults <- downloadHandler(
      
      filename = function() {
        fn<-sub( "set 1|set 2|set 3|set1|set2|set3", "", input_file()$name, ignore.case=TRUE )
        sub( "\\.csv$", "-normalized.xlsx",fn, ignore.case=TRUE )
      },
      
      content = function(file) {
        
        d_df<-req(plateData())
        d_df$UL<-NULL
        d_df$HYC01<-NULL
        d_df$HYPCON<-NULL
        d_df$UL.x<-NULL
        d_df$HYC01.x<-NULL
        d_df$HYPCON.x<-NULL
        d_df$UL.y<-NULL
        d_df$HYC01.y<-NULL
        d_df$HYPCON.y<-NULL
        # add probe names as well:
        load("Probelist_new.Rmd")
      
        bn<-look.up( colnames(d_df)[2:length(d_df)], probelist_new, to="Probe.annotation", leave.unfound=FALSE ) 
        bf<-rbind(c("Sample",bn),d_df)
        write_xlsx( bf, file)
      }
      
    )
    output$downloadResults2 <- downloadHandler(
      
      filename = function() {
        fn<-sub( "set 1|set 2|set 3|set1|set2|set3", "", input_file()$name, ignore.case=TRUE )
        sub( "\\.csv$", "-set1-normalized.xlsx",fn, ignore.case=TRUE )
      },
      
      content = function(file) {
        
        bf<-req(plateDataprobelist_set1())
        write_xlsx( bf, file)
      }
      
    )
    output$downloadResults3 <- downloadHandler(
      
      filename = function() {
        fn<-sub( "set 1|set 2|set 3|set1|set2|set3", "", input_file()$name, ignore.case=TRUE )
        sub( "\\.csv$", "-set2-normalized.xlsx",fn, ignore.case=TRUE )
      },
      
      content = function(file) {
        
        bf<-req(plateDataprobelist_set2())
        bf<-rbind(c("Sample",bn),d_df)
        write_xlsx( bf, file)
      }
      
    ) 
    output$downloadResults4 <- downloadHandler(
      
      filename = function() {
        fn<-sub( "set 1|set 2|set 3|set1|set2|set3", "", input_file()$name, ignore.case=TRUE )
        sub( "\\.csv$", "-set3-normalized.xlsx",fn, ignore.case=TRUE )
      },
      
      content = function(file) {
        
        
        bf<-req(plateDataprobelist_set3())
        write_xlsx( bf, file)
      }
      
    )
    output$excel_preview <- renderDT({
      req(interpolationToCopyNumber())
      
      datatable(
        interpolationToCopyNumber(),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    })
    output$product1 <- renderDT({
      req(plateDataprobelist_set1())
      
      datatable(
        plateDataprobelist_set1(),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    })
    output$product2 <- renderDT({
      req(plateDataprobelist_set2())
      
      datatable(
        plateDataprobelist_set2(),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    })
    output$product3 <- renderDT({
      req(plateDataprobelist_set3())
      
      datatable(
        plateDataprobelist_set3(),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    })
    output$downloadResults5 <- downloadHandler(
      
      filename = function() {
        fn<-sub( "set 1|set 2|set 3|set1|set2|set3", "", input_file()$name, ignore.case=TRUE )
        sub( "\\.csv$", "-copynumbers.xlsx",fn, ignore.case=TRUE )
      },
      
      content = function(file) {
        
        bfreshaped<-interpolationToCopyNumber()
        write_xlsx( bfreshaped, file)
      }
      
    ) 
#    output$upload_status <- renderUI({
 #     req(input$bc_file)
      
  #    tagList(
  #      lapply(seq_len(nrow(input$bc_file)), function(i) {
  #        tags$div(
  #          style = paste(
  #            "color: #198754;",
  #            "font-weight: bold;",
  #            "margin: 5px 0;"
  #          ),
  #          tags$span(
  #            style = "font-size: 18px; margin-right: 8px;",
  #            HTML("&#10003;")
  #          ),
  #          input$bc_file$name[i]
  #        )
  #      })
  #    )
  #  })
    observeEvent( input_file(), {
      if( is.null(input_file()) ) {
        shinyjs::hideElement( id="preview-button-column" )
        shinyjs::hideElement( id="preview-button-column2" )
        shinyjs::hideElement( id="preview-button-column3" )
        shinyjs::hideElement( id="preview-button-column4" )
        shinyjs::hideElement( id="download-button-column" )
        shinyjs::hideElement( id="download-button-column2" )
        shinyjs::hideElement( id="download-button-column3" )
        shinyjs::hideElement( id="download-button-column4" )
        shinyjs::hideElement( id="download-button-column5" )
      } else {
        shinyjs::showElement( id="preview-button-column" )
        shinyjs::showElement( id="preview-button-column2" )
        shinyjs::showElement( id="preview-button-column3" )
        shinyjs::showElement( id="preview-button-column4" )
        shinyjs::showElement( id="download-button-column" )
        shinyjs::showElement( id="download-button-column2" )
        shinyjs::showElement( id="download-button-column3" )
        shinyjs::showElement( id="download-button-column4" )
        shinyjs::showElement( id="download-button-column5" )
        
      }
    })  
    # get icons with eye symbols: 
    preview_sections <- list(
      
      toggle_excel_preview1 = "excel_preview_section1",
      toggle_excel_preview2 = "excel_preview_section2",
      toggle_excel_preview3 = "excel_preview_section3"
    )
    preview_visible <- reactiveValues(
      toggle_excel_preview1 = FALSE,
      toggle_excel_preview2 = FALSE,
      toggle_excel_preview3 = FALSE
    )
    
    observeEvent(plateDataprobelist_set1, {
      df <- plateDataprobelist_set1()
      req(!is.null(df), NROW(df) > 0)
      
      shinyjs::show("preview_toggle_container")
    }, ignoreInit = TRUE)
    
    ##### connection to data_processing.html:
    observeEvent(input$bc_file, {
      req(input$bc_file)
      session$sendCustomMessage("upload_done", TRUE)
    })
    observeEvent(input$resetBtn, {
      # input$bc_file is already set to NULL by Shiny.setInputValue('bc_file', null)
      # Reset any custom reactiveValues() you use here if needed
    })
    observeEvent(input$processBtn, {
      req(input$bc_file)
      session$sendCustomMessage("show_workspace", TRUE)
    })
    
    # List the uploaded file names in column(3, ...)
    output$uploaded_files_table <- renderTable({
      req(input$bc_file)
      data.frame(File = input$bc_file$name)
    })
      # Render each uploaded file name with a green checkmark icon
    # Render each uploaded file name with a green checkmark icon
    output$uploaded_files_list <- renderUI({
      req(input$bc_file)
      
      tags$ul(
        style = "list-style: none; padding-left: 0; margin-top: 10px;",
        lapply(input$bc_file$name, function(fname) {
          tags$li(
            style = "display: flex; align-items: center; gap: 8px; padding: 6px 0; border-bottom: 1px solid #f0f2f5; font-size: 14px;",
            tags$span(
              style = "display: inline-flex; align-items: center; justify-content: center; width: 20px; height: 20px; border-radius: 50%; background: #dcfce7; border: 1.5px solid #16a34a; color: #16a34a; font-weight: bold; font-size: 12px; flex-shrink: 0;",
              HTML("&#10003;")
            ),
            tags$span(fname)
          )
        })
      )
    })
    
    # Ensure Shiny renders the file list immediately even while #shiny_workspace is hidden
    outputOptions(output, "uploaded_files_list", suspendWhenHidden = FALSE)
    

  # Ensure Shiny renders the file list immediately even while #shiny_workspace is hidden
  outputOptions(output, "uploaded_files_list", suspendWhenHidden = FALSE)

    observeEvent(input$toggle_excel_preview1, {
      preview_visible(!preview_visible())
      
      if (preview_visible()) {
        shinyjs::show("excel_preview_section1")
        shinyjs::runjs(
          "
          $('#toggle_excel_preview1')
           .addClass('preview-active');
          $('#toggle_excel_preview1 i')
        .removeClass('fa-eye-slash')
        .addClass('fa-eye');
          $('#toggle_excel_preview1 i')
          .css({
           'color': '#1D9685',
           'font-weight': 'bold'
          });"
        )
      } else {
        shinyjs::hide("excel_preview_section1")
        shinyjs::runjs(
          "
          $('#toggle_excel_preview1')
           .removeClass('preview-active');
          $('#toggle_excel_preview1 i')
        .removeClass('fa-eye')
        .addClass('fa-eye-slash');
          $('#toggle_excel_preview1 i')
          .css({
           'color': '',
           'font-weight': ''
          });"
        )
      }
    })
    
    # Product 2
    observeEvent(plateDataprobelist_set2, {
      df <- plateDataprobelist_set2()
      req(!is.null(df), NROW(df) > 0)
      
      shinyjs::show("preview_toggle_container")
    }, ignoreInit = TRUE)
    observeEvent(input$toggle_excel_preview2, {
      preview_visible(!preview_visible())
      
      if (preview_visible()) {
        shinyjs::show("excel_preview_section2")
        shinyjs::runjs(
          "
          $('#toggle_excel_preview2')
           .addClass('preview-active');
          $('#toggle_excel_preview2 i')
        .removeClass('fa-eye-slash')
        .addClass('fa-eye');
          $('#toggle_excel_preview2 i')
          .css({
           'color': '#198754',
           'font-weight': 'bold'
          });"
        )
      } else {
        shinyjs::hide("excel_preview_section2")
        shinyjs::runjs(
          "
          $('#toggle_excel_preview2')
           .removeClass('preview-active');
          $('#toggle_excel_preview2 i')
        .removeClass('fa-eye')
        .addClass('fa-eye-slash');
          $('#toggle_excel_preview2 i')
          .css({
           'color': '',
           'font-weight': ''
          });"
        )
      }
    })
    # Product 3
    observeEvent(plateDataprobelist_set3, {
      df <- plateDataprobelist_set3()
      req(!is.null(df), NROW(df) > 0)
      
      shinyjs::show("preview_toggle_container")
    }, ignoreInit = TRUE)
    observeEvent(input$toggle_excel_preview3, {
      preview_visible(!preview_visible())
      
      if (preview_visible()) {
        shinyjs::show("excel_preview_section3")
        shinyjs::runjs(
          "
          $('#toggle_excel_preview3')
           .addClass('preview-active');
          $('#toggle_excel_preview3 i')
        .removeClass('fa-eye-slash')
        .addClass('fa-eye');
          $('#toggle_excel_preview3 i')
          .css({
           'color': '#198754',
           'font-weight': 'bold'
          });"
        )
      } else {
        shinyjs::hide("excel_preview_section3")
        shinyjs::runjs(
          "$('#toggle_excel_preview3')
           .removeClass('preview-active');
          $('#toggle_excel_preview3 i')
        .removeClass('fa-eye')
        .addClass('fa-eye-slash');
          $('#toggle_excel_preview3 i')
          .css({
           'color': '',
           'font-weight': ''
          });"
        )
      }
    })
    observeEvent(input$back_to_upload, {
      session$sendCustomMessage("show_landing", TRUE)
    })
    observeEvent(interpolationToCopyNumber(), {
      df <- interpolationToCopyNumber()
      req(!is.null(df), NROW(df) > 0)
      
      shinyjs::show("preview_toggle_container")
    }, ignoreInit = TRUE)
    preview_visible <- reactiveVal(FALSE)
    
    
    observeEvent(input$toggle_excel_preview, {
      preview_visible(!preview_visible())
      
      if (preview_visible()) {
        shinyjs::show("excel_preview_section")
        shinyjs::runjs("
          $('#toggle_excel_preview')
           .addClass('preview-active');
          $('#toggle_excel_preview i')
           .removeClass('fa-eye-slash')
           .addClass('fa-eye');
          $('#toggle_excel_preview i')
          .css({
           'color': '#1D9685',
           'font-weight': 'bold'
          });"
        )
      } else {
        shinyjs::hide("excel_preview_section")
        shinyjs::runjs(
          
          "$('#toggle_excel_preview')
           .removeClass('preview-active');
        $('#toggle_excel_preview i')
        .removeClass('fa-eye')
        .addClass('fa-eye-slash');
          $('#toggle_excel_preview i')
         .css({
           'color': '',
           'font-weight': ''
         });"
        )
      }
    })
   # observeEvent(input$qc_overview, {
      
  #    showModal(
  #      modalDialog(
  #        title = "QC Overview",
  #        size = "l",
  #        easyClose = TRUE,
  #        footer = modalButton("Close"),
          
   #       uiOutput("qcProductsAll"),
  #        dataTableOutput("qcTableST", width = "100%"),
  #        dataTableOutput("qcTableQCC29", width = "100%"),
  #        dataTableOutput("qcTableQCC33", width = "100%"),
  #        dataTableOutput("qcTableZymo", width = "100%"),
          
   #       tags$script(HTML("
  #      setTimeout(function() {
  #        $('.modal-dialog').last().addClass('qc-wide-modal');
  #      }, 100);
  #    "))
  #      )
  #    )
      
  #  })
    observeEvent(input$qc_overview, {
      updateCheckboxInput(
        session,
        "show_qc_page",
        value = TRUE
      )
    })
    
    observeEvent(input$qc_back, {
      updateCheckboxInput(
        session,
        "show_qc_page",
        value = FALSE
      )
    })
    observeEvent(input$pc_plot, {
      updateCheckboxInput(session, "show_pc_page", value = TRUE)
    })
    
    observeEvent(input$pc_back, {
      updateCheckboxInput(session, "show_pc_page", value = FALSE)
    })
    observeEvent(input$qc_to_pc, {
      updateCheckboxInput(session, "show_qc_page", value = FALSE)
      updateCheckboxInput(session, "show_pc_page", value = TRUE)
    })
    
    # Switch from PC Plot to QC Overview
    observeEvent(input$pc_to_qc, {
      updateCheckboxInput(session, "show_pc_page", value = FALSE)
      updateCheckboxInput(session, "show_qc_page", value = TRUE)
    })
    
   # PAGE 2
    ## Page 2: QC tables
    
  
    output$qcTableST <- DT::renderDT({
      
      df <- qctable("ST")
      cn <- names(df)
      
      # Identify product columns
      cols_1 <- cn[grepl("_1$", cn)]
      cols_2 <- cn[grepl("_2$", cn)]
      cols_3 <- cn[grepl("_3$", cn)]
      
      
      # Reorder columns by product
      df <- df[, c("Sample",cols_1, cols_2, cols_3,"Result"), drop = FALSE]
      
      # Group sizes
      n1 <- length(cols_1)
      n2 <- length(cols_2)
      n3 <- length(cols_3)
      
      
      # Build a two-row header
      group_names <- c(
        if (n1 > 0) "Product 1",
        if (n2 > 0) "Product 2",
        if (n3 > 0) "Product 3"
      )
      
      group_widths <- c(n1, n2, n3)
      group_widths <- group_widths[group_widths > 0]
      # make an empty header
      header <- htmltools::withTags(
        table(
          class = "display"
        )
      )
      
      # Find result columns that actually exist
      result_cols <- intersect(
        c("Result", "Result_1", "Result_2", "Result_3"),
        names(df)
      )

    dt <- DT::datatable(
        df,
        container = header,
        escape = FALSE, # Required to let the SVG markup render visually
        rownames = FALSE,
        options = list(
          pageLength = 10,
          # drawCallback runs automatically every time the table renders, sorts, or changes pages
          drawCallback = JS(
            "function(settings) {",
            "  $(this).find('td').each(function() {",
            "    var cell = $(this);",
            "    var txt = cell.text().trim();",
            "    var failSvg = '<svg xmlns=\"http://w3.org\" width=\"14\" height=\"14\" fill=\"currentColor\" viewBox=\"0 0 16 16\" style=\"margin-right: 6px; display: inline-block; vertical-align: -0.15em;\"><path d=\"M16 8A8 8 0 1 1 0 8a8 8 0 0 1 16 0zM5.354 4.646a.5.5 0 1 0-.708.708L7.293 8l-2.647 2.646a.5.5 0 0 0 .708.708L8 8.707l2.646 2.647a.5.5 0 0 0 .708-.708L8.707 8l2.647-2.646a.5.5 0 0 0-.708-.708L8 7.293 5.354 4.646z\"/></svg>';",
            "    var passSvg = '<svg xmlns=\"http://w3.org\" width=\"14\" height=\"14\" fill=\"currentColor\" viewBox=\"0 0 16 16\" style=\"margin-right: 6px; display: inline-block; vertical-align: -0.15em;\"><path d=\"M16 8A8 8 0 1 1 0 8a8 8 0 0 1 16 0zm-3.97-3.03a.75.75 0 0 0-1.08.022L7.477 9.417 5.384 7.323a.75.75 0 0 0-1.06 1.06L6.97 11.03a.75.75 0 0 0 1.079-.02l3.992-4.99a.75.75 0 0 0-.01-1.05z\"/></svg>';",
            "    ",
            "    // If it hasn't already been formatted, inject the icon",
            "    if (cell.find('svg').length === 0) {",
            "      if (txt === 'Fail') {",
            "        cell.html(failSvg + txt);",
            "      } else if (txt === 'Pass') {",
            "        cell.html(passSvg + txt);",
            "      }",
            "    }",
            "  });",
            "}"
          )
        )
      )
      
      # Apply colors to each result column
      for (col in result_cols) {
        dt <- dt %>%
          DT::formatStyle(
            columns = col,
            target = "cell",
            color = DT::styleEqual(
              c("Fail", "Pass"),
              c("red", "green")
            )
          )
      }
      # Color Pass and Fail cells
      for (col in result_cols) {
        dt <- dt %>%
          DT::formatStyle(
            columns = col,
            target = "cell",
            color = DT::styleEqual(
              c("Fail", "Pass"),
              c("red", "green")
            )
          )
      }
      
      # Vertical separator before Product 2
      if (n2 > 0) {
        dt <- dt %>%
          DT::formatStyle(
            columns = cols_2[1],
            target = "cell",
            borderLeft = "2px solid black"
          )
      }
      
      # Vertical separator before Product 3
      if (n3 > 0) {
        dt <- dt %>%
          DT::formatStyle(
            columns = cols_3[1],
            target = "cell",
            borderLeft = "2px solid black"
          )
      }
      
      # Vertical separator before Other columns
      if ((n1 + n2 + n3) > 0) {
        dt <- dt %>%
          DT::formatStyle(
            columns = "Result",
            target = "cell",
            borderLeft = "2px solid black"
          )
      }
      
      dt
    })
   
    output$qcTableQCC29 <- DT::renderDT({
      
      df <- qctable("QCC29")
      cn <- names(df)
      
      # Identify product columns
      cols_1 <- cn[grepl("_1$", cn)]
      cols_2 <- cn[grepl("_2$", cn)]
      cols_3 <- cn[grepl("_3$", cn)]
      
      
      # Reorder columns by product
      df <- df[, c("Sample",cols_1, cols_2, cols_3,"Result"), drop = FALSE]
      
      # Group sizes
      n1 <- length(cols_1)
      n2 <- length(cols_2)
      n3 <- length(cols_3)
      
      
      # Build a two-row header
      group_names <- c(
        if (n1 > 0) "Product 1",
        if (n2 > 0) "Product 2",
        if (n3 > 0) "Product 3"
      )
      
      group_widths <- c(n1, n2, n3)
      group_widths <- group_widths[group_widths > 0]
      # make an empty header
      header <- htmltools::withTags(
        table(
          class = "display"
        )
      )
      
      # Find result columns that actually exist
      result_cols <- intersect(
        c("Result", "Result_1", "Result_2", "Result_3"),
        names(df)
      )
      
      dt <- DT::datatable(
        df,
        container = header,
        escape = FALSE, # Required to let the SVG markup render visually
        rownames = FALSE,
        options = list(
          pageLength = 10,
          # drawCallback runs automatically every time the table renders, sorts, or changes pages
          drawCallback = JS(
            "function(settings) {",
            "  $(this).find('td').each(function() {",
            "    var cell = $(this);",
            "    var txt = cell.text().trim();",
            "    var failSvg = '<svg xmlns=\"http://w3.org\" width=\"14\" height=\"14\" fill=\"currentColor\" viewBox=\"0 0 16 16\" style=\"margin-right: 6px; display: inline-block; vertical-align: -0.15em;\"><path d=\"M16 8A8 8 0 1 1 0 8a8 8 0 0 1 16 0zM5.354 4.646a.5.5 0 1 0-.708.708L7.293 8l-2.647 2.646a.5.5 0 0 0 .708.708L8 8.707l2.646 2.647a.5.5 0 0 0 .708-.708L8.707 8l2.647-2.646a.5.5 0 0 0-.708-.708L8 7.293 5.354 4.646z\"/></svg>';",
            "    var passSvg = '<svg xmlns=\"http://w3.org\" width=\"14\" height=\"14\" fill=\"currentColor\" viewBox=\"0 0 16 16\" style=\"margin-right: 6px; display: inline-block; vertical-align: -0.15em;\"><path d=\"M16 8A8 8 0 1 1 0 8a8 8 0 0 1 16 0zm-3.97-3.03a.75.75 0 0 0-1.08.022L7.477 9.417 5.384 7.323a.75.75 0 0 0-1.06 1.06L6.97 11.03a.75.75 0 0 0 1.079-.02l3.992-4.99a.75.75 0 0 0-.01-1.05z\"/></svg>';",
            "    ",
            "    // If it hasn't already been formatted, inject the icon",
            "    if (cell.find('svg').length === 0) {",
            "      if (txt === 'Fail') {",
            "        cell.html(failSvg + txt);",
            "      } else if (txt === 'Pass') {",
            "        cell.html(passSvg + txt);",
            "      }",
            "    }",
            "  });",
            "}"
          )
        )
      )
      
      # Apply colors to each result column
      for (col in result_cols) {
        dt <- dt %>%
          DT::formatStyle(
            columns = col,
            target = "cell",
            color = DT::styleEqual(
              c("Fail", "Pass"),
              c("red", "green")
            )
          )
      }
      # Color Pass and Fail cells
      for (col in result_cols) {
        dt <- dt %>%
          DT::formatStyle(
            columns = col,
            target = "cell",
            color = DT::styleEqual(
              c("Fail", "Pass"),
              c("red", "green")
            )
          )
      }
      
      # Vertical separator before Product 2
      if (n2 > 0) {
        dt <- dt %>%
          DT::formatStyle(
            columns = cols_2[1],
            target = "cell",
            borderLeft = "2px solid black"
          )
      }
      
      # Vertical separator before Product 3
      if (n3 > 0) {
        dt <- dt %>%
          DT::formatStyle(
            columns = cols_3[1],
            target = "cell",
            borderLeft = "2px solid black"
          )
      }
      
      # Vertical separator before Other columns
      if ((n1 + n2 + n3) > 0) {
        dt <- dt %>%
          DT::formatStyle(
            columns = "Result",
            target = "cell",
            borderLeft = "2px solid black"
          )
      }
      
      dt
    })
    output$qcTableQCC33 <- DT::renderDT({
      
      df <- qctable("QCC33")
      cn <- names(df)
      
      # Identify product columns
      cols_1 <- cn[grepl("_1$", cn)]
      cols_2 <- cn[grepl("_2$", cn)]
      cols_3 <- cn[grepl("_3$", cn)]
      
      
      # Reorder columns by product
      df <- df[, c("Sample",cols_1, cols_2, cols_3,"Result"), drop = FALSE]
      
      # Group sizes
      n1 <- length(cols_1)
      n2 <- length(cols_2)
      n3 <- length(cols_3)
      
      
      # Build a two-row header
      group_names <- c(
        if (n1 > 0) "Product 1",
        if (n2 > 0) "Product 2",
        if (n3 > 0) "Product 3"
      )
      
      group_widths <- c(n1, n2, n3)
      group_widths <- group_widths[group_widths > 0]
      # make an empty header
      header <- htmltools::withTags(
        table(
          class = "display"
        )
      )
      
      # Find result columns that actually exist
      result_cols <- intersect(
        c("Result", "Result_1", "Result_2", "Result_3"),
        names(df)
      )
      
      dt <- dt <- DT::datatable(
        df,
        container = header,
        escape = FALSE, # Required to let the SVG markup render visually
        rownames = FALSE,
        options = list(
          pageLength = 10,
          # drawCallback runs automatically every time the table renders, sorts, or changes pages
          drawCallback = JS(
            "function(settings) {",
            "  $(this).find('td').each(function() {",
            "    var cell = $(this);",
            "    var txt = cell.text().trim();",
            "    var failSvg = '<svg xmlns=\"http://w3.org\" width=\"14\" height=\"14\" fill=\"currentColor\" viewBox=\"0 0 16 16\" style=\"margin-right: 6px; display: inline-block; vertical-align: -0.15em;\"><path d=\"M16 8A8 8 0 1 1 0 8a8 8 0 0 1 16 0zM5.354 4.646a.5.5 0 1 0-.708.708L7.293 8l-2.647 2.646a.5.5 0 0 0 .708.708L8 8.707l2.646 2.647a.5.5 0 0 0 .708-.708L8.707 8l2.647-2.646a.5.5 0 0 0-.708-.708L8 7.293 5.354 4.646z\"/></svg>';",
            "    var passSvg = '<svg xmlns=\"http://w3.org\" width=\"14\" height=\"14\" fill=\"currentColor\" viewBox=\"0 0 16 16\" style=\"margin-right: 6px; display: inline-block; vertical-align: -0.15em;\"><path d=\"M16 8A8 8 0 1 1 0 8a8 8 0 0 1 16 0zm-3.97-3.03a.75.75 0 0 0-1.08.022L7.477 9.417 5.384 7.323a.75.75 0 0 0-1.06 1.06L6.97 11.03a.75.75 0 0 0 1.079-.02l3.992-4.99a.75.75 0 0 0-.01-1.05z\"/></svg>';",
            "    ",
            "    // If it hasn't already been formatted, inject the icon",
            "    if (cell.find('svg').length === 0) {",
            "      if (txt === 'Fail') {",
            "        cell.html(failSvg + txt);",
            "      } else if (txt === 'Pass') {",
            "        cell.html(passSvg + txt);",
            "      }",
            "    }",
            "  });",
            "}"
          )
        )
      )
      
      # Apply colors to each result column
      for (col in result_cols) {
        dt <- dt %>%
          DT::formatStyle(
            columns = col,
            target = "cell",
            color = DT::styleEqual(
              c("Fail", "Pass"),
              c("red", "green")
            )
          )
      }
      # Color Pass and Fail cells
      for (col in result_cols) {
        dt <- dt %>%
          DT::formatStyle(
            columns = col,
            target = "cell",
            color = DT::styleEqual(
              c("Fail", "Pass"),
              c("red", "green")
            )
          )
      }
      
      # Vertical separator before Product 2
      if (n2 > 0) {
        dt <- dt %>%
          DT::formatStyle(
            columns = cols_2[1],
            target = "cell",
            borderLeft = "2px solid black"
          )
      }
      
      # Vertical separator before Product 3
      if (n3 > 0) {
        dt <- dt %>%
          DT::formatStyle(
            columns = cols_3[1],
            target = "cell",
            borderLeft = "2px solid black"
          )
      }
      
      # Vertical separator before Other columns
      if ((n1 + n2 + n3) > 0) {
        dt <- dt %>%
          DT::formatStyle(
            columns = "Result",
            target = "cell",
            borderLeft = "2px solid black"
          )
      }
      
      dt
    })
    output$qcProductsAll <- renderUI({
      
      df <- qctable("QCC29")
      cn <- names(df)
      
      # Identify product columns
      cols_1 <- cn[grepl("_1$", cn)]
      cols_2 <- cn[grepl("_2$", cn)]
      cols_3 <- cn[grepl("_3$", cn)]
      
      # Include Sample and Result columns in the layout
      n_sample <- sum(cn == "Sample")
      n_result <- sum(cn == "Result")
      
      # Count columns, including Sample and Result
      n1 <- length(cols_1)
      n2 <- length(cols_2)
      n3 <- length(cols_3)
      total <- n_sample + n1 + n2 + n3 + n_result
      # Calculate positions only for products that exist
      start <- n_sample
      headings <- list()
      
      if (n1 > 0) {
        headings[[length(headings) + 1]] <- list(
          name = "Product 1",
          center = start + n1 / 2
        )
        start <- start + n1
      }
      
      if (n2 > 0) {
        headings[[length(headings) + 1]] <- list(
          name = "Product 2",
          center = start + n2 / 2
        )
        start <- start + n2
      }
      
      if (n3 > 0) {
        headings[[length(headings) + 1]] <- list(
          name = "Product 3",
          center = start + n3 / 2
        )
        start <- start + n3
      }
      
      # Display headings
      tags$div(
        style = "position: relative; width: 100%; height: 45px;",
        lapply(headings, function(h) {
          
          # Offsets depend on which products are present
          present <- vapply(headings, function(x) x$name, character(1))
          
          extra <- if (h$name == "Product 1") {
            0
          } else if (h$name == "Product 2") {
            if ("Product 1" %in% present) 3 else 0
          } else {
            if ("Product 1" %in% present &&
                "Product 2" %in% present) {
              3
            } else if ("Product 1" %in% present ||
                       "Product 2" %in% present) {
              3
            } else {
              0
            }
          }
          
          tags$div(
            h$name,
            style = paste0(
              "position: absolute; ",
              "left: ", 100 * h$center / total + extra, "%; ",
              "transform: translateX(-50%); ",
              "font-weight: bold; font-size: 20px; ",
              "white-space: nowrap;"
            )
          )
        })
      )# Calculate positions only for products that exist
start <- n_sample
headings <- list()

if (n1 > 0) {
  headings[[length(headings) + 1]] <- list(
    name = "Product 1",
    center = start + n1 / 2
  )
  start <- start + n1
}

if (n2 > 0) {
  headings[[length(headings) + 1]] <- list(
    name = "Product 2",
    center = start + n2 / 2
  )
  start <- start + n2
}

if (n3 > 0) {
  headings[[length(headings) + 1]] <- list(
    name = "Product 3",
    center = start + n3 / 2
  )
  start <- start + n3
}

# Display headings
tags$div(
  style = "position: relative; width: 100%; height: 45px;",
  lapply(headings, function(h) {

    # Offsets depend on which products are present
    present <- vapply(headings, function(x) x$name, character(1))

    extra <- if (h$name == "Product 1") {
      0
    } else if (h$name == "Product 2") {
      if ("Product 1" %in% present) 5 else 0
    } else {
      if ("Product 1" %in% present &&
          "Product 2" %in% present) {
        5
      } else if ("Product 1" %in% present ||
                 "Product 2" %in% present) {
        10
      } else {
        0
      }
    }

    tags$div(
      h$name,
      style = paste0(
        "position: absolute; ",
        "left: ", 100 * h$center / total + extra, "%; ",
        "transform: translateX(-50%); ",
        "font-weight: bold; font-size: 20px; ",
        "white-space: nowrap;"
      )
    )
  })
)
      
    })
    
    #output$qcTableZymo<- renderTable(qctable("Zymo"), digits=1)
    
   # PAGE 3: PLOT PC PLOT FOR ALL SAMPLES: NORMALIZED DATA
    output$PCplot<-renderPlot({ 
     
      for(nr in 1:length(input$bc_file[,1])){
         
         data<-input$bc_file[[nr,'datapath']]
         Set1x1<-discovery(data,stop.at="file")
         
         Set1x <- as.matrix( Set1x1[ , grepl( probe.ibd.re(include.technical=TRUE), colnames(Set1x1) ) ] )
         rownames(Set1x) <- Set1x1$Sample
         #idx<-Set1x[,"UL"]<50
         i.qcc<-grepl("QCC23|QCC33|ST|Zymo|QCC29",Set1x1$Sample)
        
         
         Set1.hyc.i <-discovery(Set1x1,start.from="file",stop.at="qcc29",qc.check.UL=FALSE)
        
         rownames(Set1.hyc.i) <- Set1x1$Sample
         Set1.hyc.i<-Set1.hyc.i[!i.qcc,]
         if(dim(Set1.hyc.i)[1]==0){stop(paste("No patient samples in the file!"))}
         Set1.hyc.i[Set1.hyc.i < 5] <- 2.5
         
         if(nr==1){
           signal.data<-Set1.hyc.i
         }
         else{
           signal.data<-cbind(signal.data,Set1.hyc.i)
           }
      }
      
      
      pr.sd<-apply(signal.data,2,sd)
      i.sd0<-pr.sd==0
      signal.data<-signal.data[,!i.sd0]
      
      d<-as.data.frame(Set1x1[!i.qcc,]$Sample)
      colnames(d)<-"Sample"
      pca <- prcomp(log2(signal.data))
      t <- "log data"
      
      p=plot.pca.scores.12(pca,d,title = t)
      p
    })
    
    
    
    ## Page 4: Probe Annotations

    load("Probelist_new.Rmd")

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


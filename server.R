library(shiny)
source("helper.R")
source("global.R")

shinyServer(
  function(input, output, session){
    
     
#====================== Cluster Data ====================================        
    output$table <- renderTable({
      GeneChoice[[ input$datasets ]]
    })

  
#====================== Heat Map ==========================================   
    
  getTestgroup<-reactive({GeneListForR <-  GeneChoice[[input$datasets]]## This is the gene list and the student nees to load this.
                          AffyList<-GeneListForR$PROBE
                          subset_loggExprs<-merge(AffyList,gExprs,by.x="x",by.y="row.names")
                          subset_m<-melt(subset_loggExprs,id=c("x"))
                          Groups_loggexprs<-merge(subset_m,ALL_groups,by.x="variable",by.y="Var2")
                          })

  output$heatmap <- renderD3heatmap({
                                      testgroups<-subset(getTestgroup(), Norm_HR=="Other"|Norm_HR==input$checkpatient) ## Here the test groups are assigned from the 5 patient groups above. Can make this 2 variables called by the function. GSM number defines the patient.
                                      MeanGroups<-dcast(testgroups,x~Norm_HR, value.var="value",mean)
                                      Meancent<-merge(testgroups, MeanGroups, by.x="x",by.y="x")
                                      Meancent<-transform(Meancent, MeanCentOther=value-Other)
                                      Meancent$Group_GSM<-do.call(paste, c(Meancent[c("Norm_HR", "variable")], sep = "_"))
                                      Meancent<-merge(Meancent, GeneChoice[[input$datasets]], by.x="x",by.y="PROBE")
                                      Meancent$Gene_AffyID<-do.call(paste, c(Meancent[c("GENE.SYMBOL", "x")], sep = "_"))
                                      subsetMeancent<-subset(Meancent, Norm_HR==input$checkpatient,row.names=FALSE)
                                      cluster<-dcast(subsetMeancent, x~Group_GSM, value.var="MeanCentOther")
                                      cluster1<- cluster %>% select(-x)
                                      rownames(cluster1)<-cluster$x
                                      rbg <- maPalette(low="darkblue", high="red4", mid="grey", k=200)
                                      bicluster<-d3heatmap(as.matrix(cluster1), color=rbg, dendrogram="both", scale="row", k_row=input$n_rows, k_col=input$n_cols)
                                        })  
#======================================HeatMap2==========================================================================

      makeHeatMap <- reactive({
        testgroups<-subset(getTestgroup(), Norm_HR=="Other"|Norm_HR==input$checkpatient) ## Here the test groups are assigned from the 5 patient groups above. Can make this 2 variables called by the function. GSM number defines the patient.
        MeanGroups<-dcast(testgroups,x~Norm_HR, value.var="value",mean)
        Meancent<-merge(testgroups, MeanGroups, by.x="x",by.y="x")
        Meancent<-transform(Meancent, MeanCentOther=value-Other)
        Meancent$Group_GSM<-do.call(paste, c(Meancent[c("Norm_HR", "variable")], sep = "_"))
        Meancent<-merge(Meancent, GeneChoice[[input$datasets]], by.x="x",by.y="PROBE")
        Meancent$Gene_AffyID<-do.call(paste, c(Meancent[c("GENE.SYMBOL", "x")], sep = "_"))
        subsetMeancent<-subset(Meancent, Norm_HR==input$checkpatient,row.names=FALSE)
        cluster<-dcast(subsetMeancent, x~Group_GSM, value.var="MeanCentOther")
        cluster1<- cluster %>% select(-x)
        rownames(cluster1)<-cluster$x
        pairs.breaks <- seq(from=input$lower,to=input$upper,length.out=257)
        blue_red_bicluster<-heatmap.2(as.matrix(cluster1), col=bluered(256), dendrogram="both", 
                                       trace="none", breaks=pairs.breaks)
        
        blue_red_bicluster
  })
  
  
  output$heatmap2 <- renderPlot({
                                makeHeatMap()
                                      })  
  
  
  #================================Statistics==========================================================
  
    output$probeset<-renderUI({data2<- GeneChoice[[ input$datasets ]]
                               data2$PROBE<-as.character(data2$PROBE)
                               checkboxGroupInput("checkprobe", "Choose probeset", choices=data2[ ,2], select=data2[1:2,2])})

 
    output$checkdat<-renderPlot({testgroups<-subset(getTestgroup(), x==input$checkprobe & (Norm_HR==input$checkpatient|Norm_HR=="Other"))
                                  ggplot(data = testgroups, aes(x = x, y = value, fill = Norm_HR)) + geom_boxplot() + coord_flip()})


    output$Diagnostics <- renderPlot({
      testgroups<-subset(getTestgroup(), x==input$checkprobe & (Norm_HR==input$checkpatient|Norm_HR=="Other"))
      ANOVA_probes(testgroups)   
      })  

    output$Tukey <- renderPlot({
      testgroups<-subset(getTestgroup(), x==input$checkprobe & (Norm_HR==input$checkpatient|Norm_HR=="Other"))
      ANOVA_TukeyHSD(testgroups)   
      })
    
    output$Summary <- renderTable({
      testgroups<-subset(getTestgroup(), x==input$checkprobe & (Norm_HR==input$checkpatient|Norm_HR=="Other"))
      ANOVA_Summary(testgroups)   
      })
    
    output$pval <- renderTable({
      testgroups<-subset(getTestgroup(), x==input$checkprobe & (Norm_HR==input$checkpatient|Norm_HR=="Other"))
      ANOVA_pval(testgroups)   
      })

#====================== Power test ====================================    
    output$mPower <- renderPlot({
      args <- list()
      args$mu  <- input$mu
      args$sigma <- input$sigma
      args$n  <- input$n
      do.call(makePower, args)
    })
    
  }
)
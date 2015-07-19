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
    


  output$heatmap <- renderD3heatmap({GeneListForR <-  GeneChoice[[input$datasets]]## This is the gene list and the student nees to load this.
                                        AffyList<-GeneListForR$PROBE
                                        subset_loggExprs<-merge(AffyList,gExprs,by.x="x",by.y="row.names")
                                        subset_m<-melt(subset_loggExprs,id=c("x"))
                                        Groups_loggexprs<-merge(subset_m,ALL_groups,by.x="variable",by.y="Var2")
                                        testgroups<-subset(Groups_loggexprs, Norm_HR=="Other"|Norm_HR==input$checkpatient) ## Here the test groups are assigned from the 5 patient groups above. Can make this 2 variables called by the function. GSM number defines the patient.
                                        MeanGroups<-dcast(testgroups,x~Norm_HR, value.var="value",mean)
                                        Meancent<-merge(testgroups, MeanGroups, by.x="x",by.y="x")
                                        Meancent<-transform(Meancent, MeanCentOther=value-Other)
                                        Meancent$Group_GSM<-do.call(paste, c(Meancent[c("Norm_HR", "variable")], sep = "_"))
                                        Meancent<-merge(Meancent, GeneListForR, by.x="x",by.y="PROBE")
                                        Meancent$Gene_AffyID<-do.call(paste, c(Meancent[c("GENE.SYMBOL", "x")], sep = "_"))
                                        subsetMeancent<-subset(Meancent, Norm_HR==input$checkpatient,row.names=FALSE)
                                        cluster<-dcast(subsetMeancent, x~Group_GSM, value.var="MeanCentOther")
                                        cluster1<- cluster %>% select(-x)
                                        rownames(cluster1)<-cluster$x
                                        makeD3HeatMap(cluster1)  
                                        })  


    output$probeset<-renderUI({data2<- GeneChoice[[ input$datasets ]]
                               data2$PROBE<-as.character(data2$PROBE)
                               checkboxGroupInput("checkprobe", "Choose probeset", choices=data2[ ,2], select=data2[1:2,2])})

 
    output$checkdat<-renderPlot({GeneListForR <-  GeneChoice[[input$datasets]]## This is the gene list and the student nees to load this.
                                  AffyList<-GeneListForR$PROBE
                                  subset_loggExprs<-merge(AffyList,gExprs,by.x="x",by.y="row.names")
                                  subset_m<-melt(subset_loggExprs,id=c("x"))
                                  Groups_loggexprs<-merge(subset_m,ALL_groups,by.x="variable",by.y="Var2")
                                  testgroups<-subset(Groups_loggexprs, x==input$checkprobe & (Norm_HR==input$checkpatient|Norm_HR=="Other"))
                                  ggplot(data = testgroups, aes(x = x, y = value, fill = Norm_HR)) + geom_boxplot() + coord_flip()})


    output$Diagnostics <- renderPlot({
      GeneListForR <-  GeneChoice[[input$datasets]]## This is the gene list and the student nees to load this.
      AffyList<-GeneListForR$PROBE
      subset_loggExprs<-merge(AffyList,gExprs,by.x="x",by.y="row.names")
      subset_m<-melt(subset_loggExprs,id=c("x"))
      Groups_loggexprs<-merge(subset_m,ALL_groups,by.x="variable",by.y="Var2")
      testgroups<-subset(Groups_loggexprs, x==input$checkprobe & (Norm_HR==input$checkpatient|Norm_HR=="Other"))
      ANOVA_probes(testgroups)   
      })  

    output$Tukey <- renderPlot({
      GeneListForR <-  GeneChoice[[input$datasets]]## This is the gene list and the student nees to load this.
      AffyList<-GeneListForR$PROBE
      subset_loggExprs<-merge(AffyList,gExprs,by.x="x",by.y="row.names")
      subset_m<-melt(subset_loggExprs,id=c("x"))
      Groups_loggexprs<-merge(subset_m,ALL_groups,by.x="variable",by.y="Var2")
      testgroups<-subset(Groups_loggexprs, x==input$checkprobe & (Norm_HR==input$checkpatient|Norm_HR=="Other"))
      ANOVA_TukeyHSD(testgroups)   
      })
    output$Summary <- renderTable({
      GeneListForR <-  GeneChoice[[input$datasets]]## This is the gene list and the student nees to load this.
      AffyList<-GeneListForR$PROBE
      subset_loggExprs<-merge(AffyList,gExprs,by.x="x",by.y="row.names")
      subset_m<-melt(subset_loggExprs,id=c("x"))
      Groups_loggexprs<-merge(subset_m,ALL_groups,by.x="variable",by.y="Var2")
      testgroups<-subset(Groups_loggexprs, x==input$checkprobe & (Norm_HR==input$checkpatient|Norm_HR=="Other"))
      ANOVA_Summary(testgroups)   
      })
    output$pval <- renderTable({
      GeneListForR <-  GeneChoice[[input$datasets]]## This is the gene list and the student nees to load this.
      AffyList<-GeneListForR$PROBE
      subset_loggExprs<-merge(AffyList,gExprs,by.x="x",by.y="row.names")
      subset_m<-melt(subset_loggExprs,id=c("x"))
      Groups_loggexprs<-merge(subset_m,ALL_groups,by.x="variable",by.y="Var2")
      testgroups<-subset(Groups_loggexprs, x==input$checkprobe & (Norm_HR==input$checkpatient|Norm_HR=="Other"))
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
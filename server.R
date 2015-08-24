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
                                      cluster1<- cluster %>% dplyr::select(-x)
                                      rownames(cluster1)<-cluster$x
                                      makeD3HeatMap(cluster1)  
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
        cluster1<- cluster %>% dplyr::select(-x)
        rownames(cluster1)<-cluster$x
        pairs.breaks <- seq(from=input$lower,to=input$upper,length.out=257)
        blue_red_bicluster<-heatmap.2(as.matrix(cluster1), col=bluered(256), dendrogram="both", 
                                       trace="none", breaks=pairs.breaks)
        
        blue_red_bicluster
  })
  
  
  output$heatmap2 <- renderPlot({
                                makeHeatMap()
                                      })  
  
  output$bwplot_meancent<-renderPlot({testgroups<-subset(getTestgroup(), Norm_HR=="Other"|Norm_HR==input$checkpatient) ## Here the test groups are assigned from the 5 patient groups above. Can make this 2 variables called by the function. GSM number defines the patient.
                MeanGroups<-dcast(testgroups,x~Norm_HR, value.var="value",mean)
                Meancent<-merge(testgroups, MeanGroups, by.x="x",by.y="x")
                Meancent<-transform(Meancent, MeanCentOther=value-Other)
                Meancent$Group_GSM<-do.call(paste, c(Meancent[c("Norm_HR", "variable")], sep = "_"))
                Meancent<-merge(Meancent, GeneChoice[[input$datasets]], by.x="x",by.y="PROBE")
                Meancent$Gene_AffyID<-do.call(paste, c(Meancent[c("GENE.SYMBOL", "x")], sep = "_"))
                subsetMeancent<-subset(Meancent, Norm_HR==input$checkpatient,row.names=FALSE)
                bwplot(MeanCentOther~Norm_HR|Gene_AffyID, data=Meancent)
                                        })
  
  
#================================Statistics==========================================================
  
  
  output$probeset<-renderUI({data2<- GeneChoice[[ input$datasets ]]
    data2$PROBE<-as.character(data2$PROBE)
    checkboxGroupInput("checkprobe", "Choose probeset", choices=data2[ ,2], select=data2[ ,2])})
  
  mydataset<-reactive({
    identifyprobe<-input$checkprobe
    testgroups<-subset(getTestgroup(), Norm_HR=="Other"|Norm_HR==input$checkpatient) ## Here the test groups are assigned from the 5 patient groups above. Can make this 2 variables called by the function. GSM number defines the patient.
    testgroups_subset<-subset(testgroups, x %in% identifyprobe)## This allows for subseeting based on a list for rows
    colnames(testgroups_subset)[1]<-"GSM"
    colnames(testgroups_subset)[2]<-"variable"
    testgroups_ALL<-testgroups_subset
  })
  
  output$checkdat<-renderPlot({
    testgroups<-mydataset()
    ggplot(data = testgroups, aes(x = variable, y = value, fill = Norm_HR)) + geom_boxplot() 
  ##+ coord_flip()
                              })
  
  
  output$Diagnostics <- renderPlot({
    testgroups_ALL<-mydataset()
    ANOVA_residuals(testgroups_ALL)
                                  })
  
  ##output$Tukey <- renderPlot({
    ##testgroups_ALL<-mydataset()
    ##ANOVA_TukeyHSD(testgroups_ALL)   
  ##})
  
  output$Summary <- renderTable({
    testgroups_ALL<-mydataset()
    ANOVA_Summary(testgroups_ALL)   
  })
  
  output$pval <- renderTable({
    testgroups_ALL<-mydataset()
    ANOVA_pval(testgroups_ALL)   
  })
  

  output$cookdist<-renderPlot({
    testgroups_ALL<-mydataset()
    ANOVA_cooksd(testgroups_ALL) 
  })
  
  output$Stdresid<-renderPlot({
    testgroups_ALL<-mydataset()
    ANOVA_stdresiduals(testgroups_ALL)
  })
  
  output$detect_out<-renderPlot({
    testgroups_ALL<-mydataset()
    ANOVA_outlier(testgroups_ALL)
  })
  
  output$meanvals<-renderTable({
    testgroups_ALL<-mydataset()
    dcast(testgroups_ALL, variable~Norm_HR, value.var="value", mean)
  })
  
  output$samplesize<-renderTable({
    testgroups_ALL<-mydataset()
    dcast(testgroups_ALL, variable~Norm_HR, value.var="value", n_distinct)
  })
  
  
  
#====================Genome Browser====================================    
 
    output$Gene_select<-renderUI({data2a<- GeneChoice[[ input$datasets ]]
    data2a$GENE.SYMBOL<-as.character(data2a$GENE.SYMBOL)
    gene_count<-as.data.frame(tally(~GENE.SYMBOL, data2a))
    gene_count$Var1<-as.character(gene_count$Var1)
    radioButtons("checkGene", "Choose Gene", choices=gene_count[ ,1], select=gene_count[1,1])
                                  })
    
    AffyList2b<-reactive({geneKey<-input$checkGene
      AffyList2<-AnnotationDbi::select(hgu133plus2.db, keys = geneKey, columns = "PROBEID", keytype ="SYMBOL")
      Affy_Probeset<-AffyList2$PROBEID
      })
    
    returnBiomart<-reactive({geneText<-input$checkGene
    geneText<-as.character(geneText)
    r2 = getBM(attributes=c("ensembl_transcript_id","hgnc_symbol", "chromosome_name","start_position","end_position", "strand"),
              filters=c("hgnc_symbol"),
              values=list(geneText),mart=ensembl)})
    
    output$Biomart<-renderTable({r3<-returnBiomart()})
    
    chrloc<-reactive({r4<-returnBiomart()
    r4<-r4[1,3]
    r4<-paste("chr", r4, sep="")
                    })
    
    chr_l<-reactive({r5<-returnBiomart()
    r5<-r5[1,4]
    r5<-as.numeric(r5)                })
    
    chr_h<-reactive({r6<-returnBiomart()
    r6<-r6[1,5]
    r6<-as.numeric(r6)                })
    
    output$chooseProbeset<-renderUI({probesets<-AffyList2b()
      radioButtons("checkprobe_browser", "Choose probeset", choices=probesets)
      })
    
    
    output$slider1<-renderUI({sliderInput("chr_location1", label = h4("chromosome start"),
                min = chr_l(), max = chr_h(), step = 100, value = chr_l())})
    
    
    output$slider2<-renderUI({sliderInput("chr_location2", label = h4("chromosome end"),
                min = chr_l(), max = chr_h(), step = 100, value = chr_h())})
    
    chrom_probe<-reactive({ chr1<-chrloc()
                            chr1<-as.character(chr1)
                            subset(Probe_Seq_AffyIndices, L1==chr1&Probes==input$checkprobe_browser)
                            })
    
    output$chrom_start<-renderUI({chrom_probe_df<-chrom_probe()
                                  chrom_data<-chrom_probe_df["start"]
                                  checkboxGroupInput("checkstart", "Choose probe start position/s", choices=chrom_data[ ,1], select=chrom_data[ ,1])
                                  })
    
    output$chrom_table<-renderTable({chrom_probe()})
    
    output$chrom_plot<-renderPlot({chr<-chrloc()
    chr<-as.character(chr)
    from2=input$chr_location1
    to2=input$chr_location2
    ideoTrack <- IdeogramTrack(genome="hg18", chromosome=chr)
    gtr <- GenomeAxisTrack()
    txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene
    txTr <- GeneRegionTrack(txdb, genome="hg18", chromosome = chr)
    startpos<- as.numeric(input$checkstart)
    aTrack <- AnnotationTrack(start=startpos, width=25, chromosome = chr, strand="*", genome="hg18", name = "probe details", stacking = "squish")##id maps onto row names of data matrix RMANormProbe_948
    Gviz::plotTracks(c(ideoTrack, gtr, aTrack, txTr), showId = TRUE, showExonId=FALSE, from=from2, to=to2)
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
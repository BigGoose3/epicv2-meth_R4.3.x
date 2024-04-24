
## Dependencies ## 
#install.packages("shiny")
library(shiny)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(DMRcate)
library(yaml)
library(limma)
library(missMethyl)
library(Gviz)
library(minfi)
library(RColorBrewer)
library(matrixStats)
library(minfiData)
library(stringr)
library(purrr)
library(methylclock)
library(ggplot2)
library(gridExtra)
library(htmltools)

options(shiny.maxRequestSize=1000*1024^2) ## to get large zip files, should figure out how much I would want to process at once
pal <- brewer.pal(8, "Dark2")

ui <- fluidPage(
  titlePanel("DNA Methlyation Array Analysis (EPICv2 Compatible)"),
  sidebarLayout(
    sidebarPanel(
      fileInput("files", "Upload IDAT files and Sample Sheet (zip only)", accept = ".zip"),
      br(),
      actionButton("unzip","Unzip Files"),
      br(),
      selectInput("preproc_options", "Select Pre-Processing Method", 
                  c("Stratified Quantile","Functional","Illumina","Noob","SWAN","Raw"),
                  selected = "Stratified Quantile"),
      actionButton("process", "Process Data"),
      br(),
      verbatimTextOutput("unzip1"),
      tags$div(
        h4("Pipeline Citation"),
        p("This pipeline is based on the following work:"),
        p("Maksimovic J, Phipson B, Oshlack A. A cross-package Bioconductor workflow for analysing methylation array data. F1000Res. 2016 Jun 8;5:1281. doi: 10.12688/f1000research.8839.3. PMID: 27347385; PMCID: PMC4916993.")
      )),
    mainPanel(
      tabsetPanel(
        tabPanel("File Information",
                 tableOutput("zipped")),
        tabPanel("QC Plots", 
                 tags$h4("Quality Control"),
                 actionButton("qc", "Generate QC Plots"),
                 downloadButton("download_qc", label = "Download QC Report"),
                 plotOutput("qc_plots"),
                 plotOutput("normalization")),
        tabPanel("MDS Plots", 
                 tags$h4("MDS Plots"),
                 actionButton("mds", "Generate MDS Plots"),
                 plotOutput("mds1")),
        tabPanel("Diff. Methylated Positions / GSEA",
                 tags$h4("Generate DMPs"),
                 actionButton("dmp_var", "Read Variables"),
                 checkboxGroupInput("dmp_comp", "Select Variables to Compare"),
                 #verbatimTextOutput("selected_values"),
                 actionButton("dmps", "Process Differentially Methylated Positions"),
                 verbatimTextOutput("dmp_summary"),
                 selectInput("dmp_num_comp", "Which comparison would you like to generate data for?", choices = NULL),
                 downloadButton("download_dmp", "Download DMPs File (as csv)"),
                 tags$hr(),
                 tags$h4("Gene Set Enrichment Analysis"),
                 fluidRow(
                   column(width = 6,
                          radioButtons("gokegg", "Select GSEA Tool", choices = c("GO","KEGG","GSA"),
                                       selected = "GO", inline = TRUE),
                          tags$h6("If select GSA, must upload desired gene set file on right"),
                          numericInput("topGO","Select Number of Top Pathways to Display", value = 10, min = 1, step = 1)),
                   column(width = 3,
                          tags$h5("Upload Gene Set File Here (for GSA)"),
                          fileInput("uploadGSA", NULL),
                          p("Downloadable gene sets curated by Broad Molecular Institute found here: ", tags$a("WEHI Bioinformatics", href = "https://bioinf.wehi.edu.au/software/MSigDB/")))),
                 actionButton("GSEA_display", "Display Top Pathways"),
                 downloadButton("GSEA_download", "Download All Pathways"),
                 verbatimTextOutput("topGK")),
        tabPanel("Diff. Methylated Regions (DMRs)",
                 tags$h4("Generate DMRs"),
                 actionButton("dmr_var", "Read Variables"),
                 checkboxGroupInput("dmr_comp", "Select Variables to Compare (2 or more)"),
                 actionButton("dmr_select", "Select Variables"),
                 selectInput("dmr_num_comp", "Which comparison would you like to generate data for?", choices = NULL),
                 actionButton("dmrs", "Process Differentially Methylated Regions"),
                 downloadButton("download_dmr","Download DMRs file (as csv)"),
                 tags$hr(),
                 tags$h4("Plot DMRs"),
                 tags$h5("Graph using either DMR index or by exact chromosome position"),
                 tags$h6("Note: Program-driven use of the UCSC Genome Browser is limited to a maximum of one 
                                  hit every 15 seconds and no more than 5,000 hits per day."),
                 fluidRow(
                   column(width = 4,
                          numericInput("dmr_index", label = "DMR Index", value = 1),
                          actionButton("dmr_plot_index","Plot DMR using Index")),
                   column(width = 4,
                          numericInput("dmr_chrom", label = "Chromosome", value = NULL),
                          numericInput("dmr_start", label = "Region Start Position", value = NULL),
                          numericInput("dmr_end", label ="Region End Position", value = NULL),
                          actionButton("dmr_plot_coords", "Plot DMR using Coordinates"),)
                 ),
                 plotOutput("dmr_plot")),
        tabPanel("DNA Methylation Age",
                 tags$h4("Generate DNAmAge"),
                 actionButton("clocks_check", "Check Missingness of Probes"),
                 verbatimTextOutput("clocks_avail"),
                 numericInput("min_cpg", "Set Minimum Percentage of Missing CpGs to Run Clock", value = 0.20, min = 0, max = 1),
                 selectInput("clocks_choose","Choose Clocks to Evaluate",c("Horvath","Hannum", "Levine", "skinHorvath",
                                                                           "PedBE","TL"),
                             selected = c("Horvath", "skinHorvath","Levine","TL"), multiple = TRUE),
                 actionButton("gen_age", "Generate DNAm Age for Selected Clocks"),
                 verbatimTextOutput("DNAmAge"),
                 actionButton("gen_age_plots", "Generate DNAmAge Plots"),
                 plotOutput("age_plots"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Unzipping files on click of button and then storing them as a variable
  unzipping <- reactive({
    message("Unzipping...") # this message only pops up in the console, not on the page
    output$zipped <- renderTable({
      unzip(input$files$datapath, list = FALSE)
    })
    message("Unzipped")
  }) %>% bindEvent(input$unzip)
  
  output$unzip1 <- renderPrint({
    if (input$unzip)
      unzipping()
    
  })
  
  # set necessary reactives 
  ann_v2 <- reactive({
    getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  })
  
  ann_v2Sub <- reactive({
    req(ann_v2())
    req(mVals())
    ann_v2()[match(rownames(mVals()),ann_v2()$Name),
           c(1:4,12:19,24:ncol(ann_v2()))]
  })
  
  bVals <- reactive({
    req(mSetSq())
    getBeta(mSetSq())
  })
  
  mVals <- reactive({
    req(mSetSq())
    getM(mSetSq())
  })
  
  mSetSq <- reactive({
    req(rgSet)
    data <- NULL
    if (input$preproc_options == "Stratified Quantile"){
      data <- preprocessQuantile(rgSet)
    } else if (input$preproc_options == "Functional"){
      data <- preprocessFunnorm(rgSet)
    } else if (input$preproc_options == "Illumina"){
      data <- preprocessIllumina(rgSet)
    } else if (input$preproc_options == "Noob"){
      data <- preprocessNoob(rgSet)
    } else if (input$preproc_options == "SWAN"){
      data <- preprocessSWAN(rgSet)
    } else if (input$preproc_options == "Raw"){
      data <- preprocessRaw(rgSet)
    }
    
    return(data)
  })
  
  # Process Data 
  observe({
    message("Processing Data..")
    targets <- read.metharray.sheet(getwd(), pattern = "Sample.*Sheet.*csv")
        
    # edit incomplete sample sheets
    missing_pool <- is.na(targets$Pool_ID)
    missing_group <- is.na(targets$Sample_Group)
    targets$Pool_ID[missing_pool] <- 1
    targets$Sample_Group[missing_group] <- 1
    targets$Pool_ID <- gsub(" ", "_", targets$Pool_ID)
    targets$Sample_Group <- gsub(" ", "_", targets$Sample_Group)
    targets$id <- paste(targets$Slide, targets$Array, sep = "_")
    
    #complete processing
    rgSet <- read.metharray.exp(targets = targets)
    rgSet@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
    detP <<- detectionP(rgSet)
    keep <- colMeans(detP) < 0.05 ## remove those probes with poor quality signal (typically have high p-value)
    rgSet <<- rgSet[, keep]
    targets <<- targets[keep, ]
    message("Data Processed")
  }) %>% bindEvent(input$process)

  #### Generate QC / Normalization plots and generate pdf report ####
  observe({
    message("Generating  Plots")
    
    ## QC Plots
    output$qc_plots <- renderPlot({
      par(mfrow = c(1, 2))
      req(detP)
      barplot(colMeans(detP), col = pal[factor(targets$Sample_Group)], las = 2, 
            cex.names = 0.8, ylab = "Mean detection p-values")
      abline(h = 0.05, col = "red")
      legend("topleft", legend = levels(factor(targets$Sample_Group)), fill = pal,
           bg = "white")
      barplot(colMeans(detP), col = pal[factor(targets$Sample_Group)], las = 2, 
            cex.names = 0.8, ylim = c(0, 0.002), ylab = "Mean detection p-values")
      abline(h = 0.05, col = "red")
      legend("topleft", legend = levels(factor(targets$Sample_Group)), fill = pal, 
           bg = "white")
    })
    
    ## Normalization Visualization
    output$normalization <- renderPlot({
      pal <- brewer.pal(8, "Dark2")
      par(mfrow = c(1, 2))
      req(rgSet)
      densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
      legend("top", legend = levels(factor(targets$Sample_Group)), 
           text.col=brewer.pal(8,"Dark2"))
      densityPlot(getBeta(mSetSq()), sampGroups=targets$Sample_Group,
                main="Normalized", legend=FALSE)
      legend("top", legend = levels(factor(targets$Sample_Group)), 
           text.col=brewer.pal(8,"Dark2"))
    })
    message("Plots Generated")
    
  }) %>% bindEvent(input$qc)

  
  ## Generate Downloadable QC Report
  output$download_qc <- downloadHandler(
    filename = function(){
      paste0("qcReport",".pdf")
    },
    
    content = function(file){
      qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group,
               pdf = file)
    }
  )
  
  ## MDS plots 
  observe({
    message("Generating MDS Plots")
    
    ## MDS Plots (Two components)
    output$mds1 <- renderPlot({
      req(mSetSq())
      par(mfrow=c(1,2))
      plotMDS(getM(mSetSq()), top=1000, gene.selection="common", 
              col=pal[factor(targets$Sample_Group)])
      legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
             bg="white", cex=0.7)
      
      plotMDS(getM(mSetSq()), top=1000, gene.selection="common",  
              col=pal[factor(targets$Pool_ID)])
      legend("top", legend=levels(factor(targets$Pool_ID)), text.col=pal,
             bg="white", cex=0.7)
    })
  }) %>% bindEvent(input$mds)
  
  #### DM Overlap ####
  design <- reactive({
    req(targets)
    cellType <- factor(targets$Sample_Group)
    individual <- factor(targets$Pool_ID)
    # Check the number of levels in cellType and individual
    n_cellType <- nlevels(cellType)
    n_individual <- nlevels(individual)
    # Create design1 matrix with only factors having 2 or more levels (will ignore others)
    if (n_cellType >= 2 & n_individual >= 2) {
      design <- model.matrix(~0 + cellType + individual, data = targets)
    } else if (n_cellType >= 2) {
      design <- model.matrix(~0 + cellType, data = targets)
    } else if (n_individual >= 2) {
      design <- model.matrix(~0 + individual, data = targets)
    } else {
      design <- NULL
    }
  colnames(design) <- c(levels(cellType),levels(individual)[-1])
  
  return(design)
  })
  

  #### Differentially Methylated Positions / GSEA ####
  observe({
    req(targets)
    group_val <- unique(targets$Sample_Group)
    pool_val <- unique(targets$Pool_ID)
    compvar <- c(group_val, pool_val)
    updateCheckboxGroupInput(session,'dmp_comp',
                             choices = compvar, selected = NULL)
  }) %>% bindEvent(input$dmp_var)
  
  combos_dmp <- reactive({
    compvar_sel <- c(input$dmp_comp)
    var <- str_split(input$dmp_comp,",")
    combn(var, 2, function(x) paste(x, collapse = "-"))
  })
  
  observe({
    updateSelectInput(inputId = "dmp_num_comp", choices = combos_dmp())
  }) %>% bindEvent(input$dmps)
  
  observe({
    message("Obtaining Differentially Methylated Positions")
    
    # Prepare Comparisons to be made from Checkbox
    req(design())
    req(mVals())
    fit <- lmFit(mVals(), design())
    contMatrix <- makeContrasts(contrasts=combos_dmp(),
                               levels=design())
    message("Matrix Generated")
    
    fit2 <- contrasts.fit(fit, contMatrix)
    fit2 <<- eBayes(fit2)
    
    message("Fit Created")
    
    output$dmp_summary <- renderPrint({
      req(fit2)
      summary(decideTests(fit2))
    })
    
  }) %>% bindEvent(input$dmps)
  
  DMPs <-reactive({
    # get the table of results for whatever contrast you choose
    dmps <- topTable(fit2, num=Inf, coef=input$dmp_num_comp, genelist=ann_v2Sub())
    dmps <- as(dmps, "data.frame")
    return(dmps)
  }) 
  
  output$download_dmp <- downloadHandler(
    filename = function(){
      paste0(input$dmp_num_comp, "DMP", ".csv")
    },
    
    content = function(file){
      write.table(DMPs(), file, sep=",")
    }
  )
  
  ## DNA Gene Set Enrichment Analysis (GSEA) ##
  gst <- reactive({
    DMPs <- DMPs()
    DMPs$Name <- gsub("_BC11","",as.character(DMPs$Name))
    DMPs$Name <- gsub("_TC11","",as.character(DMPs$Name))
    DMPs$Name <- gsub("_BC21","",as.character(DMPs$Name))
    DMPs$Name <- gsub("_TC21","",as.character(DMPs$Name)) 
    sigCpGs <- DMPs$Name[DMPs$adj.P.Val<0.05]
    all <- DMPs$Name
    gst <- gometh(sig.cpg=sigCpGs, all.cpg=all, collection = input$gokegg, plot.bias=FALSE,
                  anno=ann_v2())
    message("GSEA Complete")
    colnames(gst)[colnames(gst) == "ONTOLOGY"] <- "Ont"
    colnames(gst)[colnames(gst) == "TERM"] <- "Term"
    colnames(gst)[colnames(gst) == "Description"] <- "Pathway"
    return(gst)
  })
  
  gsa_coll <- reactive({
    req(input$uploadGSA)
    load(input$uploadGSA$datapath)
    loaded <- ls()
    gsa_coll <- get(loaded[1])
    return(gsa_coll)
  })
  
  gsa <- reactive({
    req(gsa_coll)
    DMPs <- DMPs()
    DMPs$Name <- gsub("_BC11","",as.character(DMPs$Name))
    DMPs$Name <- gsub("_TC11","",as.character(DMPs$Name))
    DMPs$Name <- gsub("_BC21","",as.character(DMPs$Name))
    DMPs$Name <- gsub("_TC21","",as.character(DMPs$Name)) 
    sigCpGs <- DMPs$Name[DMPs$adj.P.Val<0.05]
    all <- DMPs$Name
    gsa <- gsameth(sig.cpg=sigCpGs, all.cpg=all, collection=gsa_coll())
    return(gsa)
  })
  
  observe({

    if (input$gokegg == "GO"){
      req(gst())
      topGK <- topGO(gst(), number=input$topGO)
      } else if (input$gokegg == "KEGG"){
        req(gst())
        topGK <- topKEGG(gst(), number = input$topGO)
      } else if (input$gokegg == "GSA"){
        req(gsa())
        topGK <- topGSA(gsa(), number = input$topGO)
      }
    
    output$topGK <- renderPrint({
      return(topGK)
    })
  }) %>% bindEvent(input$GSEA_display)
  
  
  
  output$GSEA_download <- downloadHandler(
    filename = function(){
      paste0(input$dmp_num_comp, "_GSEA_", input$gokegg, ".csv")
    },
    
    content = function(file){
    if (input$gokegg %in% c("GO", "KEGG")) {
      req(gst())
      write.table(gst(), file, sep = ",")
    } else if (input$gokegg == "GSA") {
      req(gsa())
      write.table(gsa(), file, sep = ",")
    }
  }
    )
  
  #### Differentially Methylated Regions ####
  observe({
    req(targets)
    group_val <- unique(targets$Sample_Group)
    pool_val <- unique(targets$Pool_ID)
    compvar <- c(group_val, pool_val)
    updateCheckboxGroupInput(session,'dmr_comp',
                             choices = compvar, selected = NULL)
  }) %>% bindEvent(input$dmr_var)
  
  combos_dmr <- reactive({
    compvar_sel <- c(input$dmr_comp)
    var <- str_split(input$dmr_comp,",")
    combn(var, 2, function(x) paste(x, collapse = "-"))
  })
  
  observe({
    updateSelectInput(inputId = "dmr_num_comp", choices = combos_dmr())
  }) %>% bindEvent(input$dmr_select)
  
  observe({
    message("Organizing Metadata for DMRs")

    # Prepare Comparisons to be made from Checkbox
    contMatrix <- makeContrasts(contrasts=combos_dmr(),
                                levels=design())
    message("Matrix Generated")
    
    # Generate Annotation
    message("Generating DMRs")
    myAnnotation <- cpg.annotate("array", mSetSq(), design = design(), 
                                 contrasts = TRUE, cont.matrix = contMatrix, 
                                 coef = input$dmr_num_comp, fdr = 0.05) 
    DMRs <<- dmrcate(myAnnotation, lambda = 1000, C = 2)
    results.ranges <<- extractRanges(DMRs, genome = "hg38")
    results.ranges.df <<- as(results.ranges, "data.frame")
    
    message("DMR Annotation Created")
    
  }) %>% bindEvent(input$dmrs)
  
  output$download_dmr <- downloadHandler(
    filename = function(){
      paste0(input$dmr_num_comp, "DMR", ".csv")
    },
    
    content = function(file){
      write.table(results.ranges.df, file, sep=",")
    }
  )
  
  # DMR Index
  observe({
    # indicate which genome is being used (could alter later)
    gen <- "hg38"
    # the index of the DMR that we will plot 
    dmrIndex <- input$dmr_index
    # extract chromosome number and location from DMR results 
    coords <- strsplit2(DMRs@coord[dmrIndex],":")
    chrom <- coords[1]
    start <- as.numeric(strsplit2(coords[2],"-")[1])
    end <- as.numeric(strsplit2(coords[2],"-")[2])
    
    # add 25% extra space to plot
    minbase <- start - (0.25*(end-start))
    maxbase <- end + (0.25*(end-start))
    
    message("Calculating iTrack")
    iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, height = 0.5)
    message("Calculating gTrack")
    gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
    message("Calculating rTrack")
    rTrack <- UcscTrack(genome=gen, chromosome=chrom, track="All GENCODE V44", 
                        from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                        rstarts="exonStarts", rends="exonEnds", gene="name", 
                        symbol="name2", transcript="name", strand="strand", 
                        fill="darkblue",stacking="squish", name="RefSeq", 
                        showId=TRUE, geneSymbol=TRUE)
    
    ann_v2Ord <- ann_v2Sub()[order(ann_v2Sub()$chr,ann_v2Sub()$pos),]
    
    bValsOrd <- bVals()[match(ann_v2Ord$Name,rownames(bVals())),]

    cpgData <- GRanges(seqnames=Rle(ann_v2Ord$chr),
                       ranges=IRanges(start=ann_v2Ord$pos, end=ann_v2Ord$pos),
                       strand=Rle(rep("*",nrow(ann_v2Ord))),
                       betas=bValsOrd)
    cpgData <- subsetByOverlaps(cpgData, results.ranges[dmrIndex])

    message("Calculating MethTrack")
    methTrack <- DataTrack(range=cpgData, groups=targets$Sample_Group,genome = gen,
                           chromosome=chrom, ylim=c(-0.05,1.05), col=pal,
                           type=c("a","p"), name="DNA Meth.\n(beta value)",
                           background.panel="white", legend=TRUE, cex.title=0.8,
                           cex.axis=0.8, cex.legend=0.8)
    message("Calculating dmrTrack")
    dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                                chromosome=chrom,fill="darkred")
    tracks <- list(iTrack, gTrack, methTrack, dmrTrack, rTrack)
    output$dmr_plot <- renderPlot({
      req(tracks)
      sizes <- c(2,2,5,2,2) # set up the relative sizes of the tracks
      try(plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
                 add35=TRUE, grid=TRUE, lty.grid=3, sizes=sizes, length(tracks)))

    })
    
  }) %>% bindEvent(input$dmr_plot_index)
  
  # Coordinates 
  observe({
    gen <- "hg38"
    chrom <- input$dmr_chrom
    start <- input$dmr_start
    end <- input$dmr_end
    
    # add 25% extra space to plot
    minbase <- start - (0.25*(end-start))
    maxbase <- end + (0.25*(end-start))
    
    message("Calculating iTrack")
    iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, height = 0.5)
    message("Calculating gTrack")
    gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
    message("Calculating rTrack")
    rTrack <- UcscTrack(genome=gen, chromosome=chrom, track="All GENCODE V44", 
                        from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                        rstarts="exonStarts", rends="exonEnds", gene="name", 
                        symbol="name2", transcript="name", strand="strand", 
                        fill="darkblue",stacking="squish", name="RefSeq", 
                        showId=TRUE, geneSymbol=TRUE)
    ann_v2Ord <- ann_v2Sub()[order(ann_v2Sub()$chr,ann_v2Sub()$pos),]
    
    bValsOrd <- bVals()[match(ann_v2Ord$Name,rownames(bVals())),]
    
    cpgData <- GRanges(seqnames=Rle(ann_v2Ord$chr),
                       ranges=IRanges(start=ann_v2Ord$pos, end=ann_v2Ord$pos),
                       strand=Rle(rep("*",nrow(ann_v2Ord))),
                       betas=bValsOrd)
    cpgData <- subsetByOverlaps(cpgData, results.ranges[dmrIndex])
    message("Calculating MethTrack")
    methTrack <- DataTrack(range=cpgData, groups=targets$Sample_Group,genome = gen,
                           chromosome=chrom, ylim=c(-0.05,1.05), col=pal,
                           type=c("a","p"), name="DNA Meth.\n(beta value)",
                           background.panel="white", legend=TRUE, cex.title=0.8,
                           cex.axis=0.8, cex.legend=0.8)
    message("Calculating dmrTrack")
    dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                                chromosome=chrom,fill="darkred")
    tracks <- list(iTrack, gTrack, methTrack, dmrTrack, rTrack)
    output$dmr_plot <- renderPlot({
      req(tracks)
      sizes <- c(2,2,5,2,2) # set up the relative sizes of the tracks
      plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
                 add35=TRUE, grid=TRUE, lty.grid=3, sizes=sizes, length(tracks))
      
    })
    
  }) %>% bindEvent(input$dmr_plot_coords)
  

  #### DNA Methylation Age Calculation ####
  bVals.df <- reactive({
    req(bVals())
    bVals_df <- as.data.frame(bVals())
    bVals_df <- tibble::rownames_to_column(bVals_df, "Probe") 
    # Won't run if there are no weird additional abbrevaitions on the probes
    # Still not entirely sure why these new EPICv2 has this 
    bVals_df$Probe <- gsub("_BC11","",as.character(bVals_df$Probe))
    bVals_df$Probe <- gsub("_TC11","",as.character(bVals_df$Probe))
    bVals_df$Probe <- gsub("_BC21","",as.character(bVals_df$Probe))
    bVals_df$Probe <- gsub("_TC21","",as.character(bVals_df$Probe)) 
    return(bVals_df)
  })
  
  # Set up and calculate cpg missingness
  observe({
    message("Check Clocks")
    output$clocks_avail <- renderPrint({
      cpgs.missing <- checkClocks(bVals.df())
    })
    message("Clocks Checked")
  }) %>% bindEvent(input$clocks_check)

  # Calculate DNAmAge, see if I can output that too for each in a table format?
  observe({
    age <- DNAmAge(bVals.df(), clocks = input$clocks_choose, min.perc = input$min_cpg)
    message("Age Generated")

    age <- as.data.frame(age)
    targets_sub <- targets[,c("Sample_Plate", "Sample_Group", "Pool_ID", "Sample_Name", "id")]
    # Identify columns with non-NA values
    cols_with_vars <- colnames(targets_sub)[colSums(!is.na(targets_sub)) > 0]
    targets_subset <- targets_sub[, cols_with_vars]
    # Merge 'age' dataframe with the subset of 'targets' dataframe
    age.df <<- merge(age, targets_subset, by.x = "id", by.y = "id", all.x = TRUE)
    
    output$DNAmAge <- renderPrint({
      return(age.df)
    })
  }) %>% bindEvent(input$gen_age)
  
  # Plot DNAmAge
  observe({
    message("Generating Plots")
    plot_list <- list()
    for(clock_choice in input$clocks_choose){
      plot <- ggplot(age.df, aes(x = Sample_Name, y = .data[[clock_choice]])) + 
        geom_point(size = 1.5) + 
        labs(title = clock_choice, y = "DNAm") +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 
      plot_list[[clock_choice]] <- plot
    }
    
    message("Gen2")
    output$age_plots <- renderPlot({
      req(plot_list)
      grobs <- lapply(plot_list, ggplotGrob)
      grid.arrange(grobs = grobs, ncol = 2)
    }, res = 96)
    
    message("Plots Generated")

  }) %>% bindEvent(input$gen_age_plots)

  
}  

shinyApp(ui, server)

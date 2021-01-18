  #
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# 
#
#    
#
library(shinythemes)
library(shiny)
library(DT)
library(dplyr)
library(stringr)
library(ggplot2)
library(tm)
library(wordcloud)
  
options(DT.options = list(pageLength = 5))

# function definitions ----------------------------------------------------

filtery_by_tissue <- function(edges_df, tissue)
{
  temp <- edges_df %>% filter(str_detect(Tissues, tissue))
}


go_analysis <- function(targeted_genes , go_terms_with_genes, num_background_genes, min_num_genes_per_term, go_category){
  # inputs: 
  # 1. targeted_genes: is a list of direct and indirect genes targeted by a miRNA
  # 2. go_terms_with_genes: a dataframe with GO terms/pathways and genes belonging to this term/pathway
  # 3. num_background_genes: number of total background genes
  # 4. min_num_genes_per_term: threshold for minimum number of genes per GO term
  # 5. go_category:  "biological_process" OR "molecular_function" OR "cellular_component"
  
  results <- go_terms_with_genes
  
  temp <-  c(1:length(results[[1]]))
  g <- rep("", length(results[[1]]))
  
  for (i in temp){
    
    inter <- intersect( targeted_genes , trimws(as.list(strsplit(results$genes[i], ","))[[1]]))
    temp[i] <- length(inter)
    g[i]    <-  paste((inter),collapse=", ")    
    
  }
  
  results$num_intersection  <- temp
  results$targeted_genes  <- g
  results$p.value  <-  1 - phyper(q = results$num_intersection  -1 , 
                                  m = results$num_genes,n = num_background_genes - results$num_genes , 
                                  k = length(targeted_genes))
  
  results <- results %>% filter( GO.domain == go_category )
  results <- results %>% filter(num_genes >= min_num_genes_per_term )
  
  results <- results %>% arrange(p.value)
  results$adj.p.value <- p.adjust(results$p.value , method ="BH" )
  results <- results %>% dplyr::select( - c(genes)   )
  
  
  return(results)

}




# load preprocessed data --------------------------------------------------
#read list of human transcription factors
hsa_TFs <- readRDS("data/hsa_TF.rds")

# list of tissues
tissues <- read.table("data/tissue_mapping.csv",strip.white = F,sep = ",",quote = "", header = T)

# gene_GO is df with genes and their assoicated go terms
go_terms_with_genes <- as.data.frame(readRDS('data/go_terms_with_genes.rds'))
#go_terms_with_genes$GO.domain <- as.character(go_terms_with_genes$GO.domain)
#go_terms <- go_terms_with_genes %>% filter(GO.domain == input$go_category) 

#gene_name_id is df with genes names and their Ensembl IDs (only genes with Go terms are included)
gene_name_id <- readRDS("data/gene_name_id.rds")

# read miRNA direct targets from TargetScan 7.2
miRNA_direct_targets <- readRDS("data/targetScan_human.rds")

# reading tissue-specific TF-targets genes associations from "tissue-specific gene regulatory networks"
ts_edges <- readRDS("data/tissue_specific_edges.rds")

ts_edges$Tissues  <- as.character(ts_edges$Tissues)

# add gene name to target genes in ts_edges and keep only genes with GO terms 
ts_edges <- ts_edges %>% inner_join(gene_name_id , by = c("TargetGene" = "Gene.stable.ID"))
colnames(ts_edges)[6] <- "symbol" 
colnames(ts_edges)[2] <- "target.gene" 
ts_edges$symbol  <- as.character(ts_edges$symbol)

# get number of background genes
ts_genes <- data.frame(Gene = unique(gene_name_id$Gene.name))
num_background_genes <- length(unique(ts_genes$Gene))

# app  UI -----------------------------------------------------------------
ui <- fluidPage( theme = shinytheme("cosmo") , 
                 
                 titlePanel("miRinGO:Prediction of biological processes targeted by human microRNAs" , ),
                 
                   sidebarPanel(
                     selectInput("t",
                                 "tissue type:",
                                 choices = c(levels(tissues$broad_tissue)),
                                 selected = "" 
                                 ),
                     
                     selectInput("miRNAs",
                                 "Select input miRNAs:",
                                 choices = miRNA_direct_targets$miRNA,
                                 multiple = T
                     ),
                     selectInput("go_category",
                                 "Select GO category:",
                                 choices = levels(go_terms_with_genes$GO.domain ),
                                 selected = "biological_process",
                                 multiple = F
                     ),
                     selectInput("mode",
                                 "Direct/Indirect targeting:",
                                 choices = c("direct","indirect"),
                                 selected = "indirect" 
                     ),
                     
                     numericInput("max_TF" , "Enter percentage of top-ranked targets (sorted by TargetScan context++ score)" , 
                                  100 , min = 20 , 100 , step = 20) ,
                     
                     numericInput("min_genes" , "Enter minimum number of genes per GO term/pathway" , 5) ,
                     numericInput("top_terms" , "Enter number of GO terms/pathways for visualization" ,10) ,
                     
                     actionButton("run", "Run Analysis!" ),
                     downloadButton("downloadData", "Download Results")
                   ),
                   
                   mainPanel(
                     
                    tabsetPanel(
                      tabPanel("Enrichment Analysis", 
                        fluidRow(
                          column(
                              dataTableOutput(outputId = "table" ), width = 6)
                                ) , 
                                    icon = icon("table" ,"fa-3x")
                              ) ,
                            # visualization page
                              tabPanel( "visualization",icon = icon("bar-chart-o","fa-3x"), 
                                        h2(strong("Bar plot")),
                                        plotOutput("barplot"),
                                        h2(strong("Word cloud")),
                                        plotOutput("wordcloud")
                                        
                                        )
                            )
                      )
)


server <- function(input, output) {
  
  tissue_specific_edges <- eventReactive( 
    {
      input$t

    }
    ,{ 
    
  ts_edges_filtered <- ts_edges 
  
  if (input$t != "all")
  {
    edges_filtered_by_tissue <- filtery_by_tissue(ts_edges_filtered , input$t ) 
    
  }else {
    edges_filtered_by_tissue <- ts_edges_filtered
  }
  
  ts_edges_summary <- edges_filtered_by_tissue %>% group_by(TF) %>% 
    summarise(number_targets = n(), genes = paste((symbol),collapse=", "))
  
  })
  
 
  
  # get miRNA targets for a miRNA family
  get_indirect_targets <- eventReactive( 
    { 
      input$t
      input$miRNAs
      input$max_TF
      input$mode
    } ,
            
  {
  
  direct_targets <- list()
  indirect_targets <- list()
 target_TFs_with_targets<- list()
  selected_miRNAs <- miRNA_direct_targets %>% filter(miRNA %in% input$miRNAs ) 
    
  for (i in c(1:length(selected_miRNAs$miRNA)))
  {

  target_genes <-  as.data.frame(trimws(as.list(strsplit(selected_miRNAs[i,]$genes, "," ) )[[1]]))
  names(target_genes) <- c("symbol")
  target_genes$symbol <- as.character(target_genes$symbol)
  
  ts_target_genes  <- target_genes %>%  filter( symbol %in% ts_genes$Gene) 
  
  target_TFs  <- ts_target_genes  %>%  filter( symbol %in% hsa_TFs$symbol) %>% dplyr::select( symbol)
  
  target_TFs_with_targets <- inner_join(target_TFs,tissue_specific_edges() ,by = c("symbol" = "TF"))
  
  # join all indirect targets in one list
  tt <- list()
  tt <- target_TFs_with_targets$symbol 
  tt <- tt[1:ceiling (input$max_TF * length(tt) / 100)  ]
  

  for (j in c(1:floor(length(tt)))){
    t1 <- trimws(as.list(strsplit(target_TFs_with_targets[j,]$genes, ","))[[1]])
    t1 <- intersect(t1,  ts_genes$Gene)
    tt <- union(tt ,  t1)
  }
  

  
  if (i==1) 
  {
    direct_targets <- target_genes$symbol
    indirect_targets <- tt
  }  else  {
    indirect_targets <- intersect(indirect_targets , tt)
    direct_targets <- intersect( direct_targets , target_genes$symbol)
  }
  
  }
  
  direct_targets <- direct_targets[1:ceiling (input$max_TF * length(direct_targets) / 100)  ]
  if (input$mode == "indirect")
  {
    indirect_targets
  } else {
    if (input$mode == "direct" )
    {
      direct_targets
    } else {
      
      both_targets <- union( indirect_targets , direct_targets)
    }
    
    
  }
  
  
  })

  enrichment <- eventReactive( input$run , { 
        

        results <- go_analysis(get_indirect_targets() ,
                               go_terms_with_genes , num_background_genes , input$min_genes , input$go_category) 
    
    })
  
  
  
  output$table <- DT::renderDataTable({ 
    x <- enrichment()
    x <- x %>% dplyr::select( - c(targeted_genes)   )

    #add hyperlinks to QuickGO GO term page
    quickGO <- "https://www.ebi.ac.uk/QuickGO/term/"
    x <- x %>% mutate( GO.term.accession =  paste0("<a href=",quickGO,GO.term.accession,
                                                   " target='_blank'>",GO.term.accession,"</a>") )
    
    } , rownames =F , escape = FALSE)
  
  # Top-K GOterms/pathways to be visualized 
    top_k <- reactive( { 
      
        res <-  enrichment()
        res <- res[1:input$top_terms,]
    })
    
  output$barplot <- renderPlot(

      
    ggplot(top_k(), aes(y = -log10(adj.p.value) , x = reorder(GO.term.name , -adj.p.value) , 
                  fill =num_intersection / num_genes ) ) +
      theme_bw() +
      theme(axis.text = element_text(size = 15) , 
            axis.title = element_text(size = 20), 
            axis.title.y = element_blank() ) + 
      geom_bar(stat = "identity" ) +
      geom_hline(yintercept = -log10(0.05) ,linetype = 'dashed' , size = 1)+
      scale_fill_gradient(low="red", high="blue") +
      coord_flip() 
    
  )  
  
  wc <- reactive(
    {
      # get enriched GO terms
      terms <- as.character(top_k()$GO.term.name)
      
      # remove uneccessary words
      myCorpus = Corpus(VectorSource(terms))
      myCorpus = tm_map(myCorpus, content_transformer(tolower))
      myCorpus = tm_map(myCorpus, removePunctuation)
      myCorpus = tm_map(myCorpus, removeNumbers)
      myCorpus = tm_map(myCorpus, removeWords,
                        c("process", "cell","regulation", "pathway","rna","dna","polymerase","via","protein",
                          "positive", "negative", "system", "transcription","mirna"))
      
      myDTM = TermDocumentMatrix(myCorpus,
                                 control = list(minWordLength = 1))
      
      m = as.matrix(myDTM)
      
      sort(rowSums(m), decreasing = TRUE)
    }
  )
  
  output$wordcloud <- renderPlot(
    
    wordcloud(names(wc()), wc(), scale=c(4,0.3),
                  min.freq = 2, max.words=50,
                  colors=brewer.pal(8, "Dark2"))
  )
  
  output$downloadData <- downloadHandler({
   filename = function() {
    paste0(c(input$miRNAs ,input$mode ,  "go_terms.txt") , collapse = "_" )
   }
  },
  content = function(filename) {
     write.table(enrichment(), filename , row.names = FALSE , quote = F , sep = '\t')
   }
  )
}


# Run the application 
shinyApp(ui = ui, server = server)



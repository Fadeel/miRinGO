library(ggplot2)
library(tm)
library(wordcloud)
library(dplyr)

d <- read.table("C:\\Users\\fadil\\Downloads\\1.txt",sep = "\t" , as.is = T , header = TRUE )

x <- d %>% slice(100:150) %>% select(GO.term.name)
# wordcloud

terms <- as.character(x$GO.term.name)

myCorpus = Corpus(VectorSource(terms))
myCorpus = tm_map(myCorpus, content_transformer(tolower))
myCorpus = tm_map(myCorpus, removePunctuation)
myCorpus = tm_map(myCorpus, removeNumbers)
myCorpus = tm_map(myCorpus, removeWords,
                  c(stopwords("SMART"), "process", "cell","regulation", "pathway",
                    "positive", "negative", "system", "transcription","miRNA"))

myDTM = TermDocumentMatrix(myCorpus,
                           control = list(minWordLength = 3))

m = as.matrix(myDTM)

v = sort(rowSums(m), decreasing = TRUE)

wordcloud_rep <- repeatable(wordcloud)

wordcloud_rep(names(v), v, scale=c(4,0.3),
              min.freq = 2, max.words=50,
              colors=brewer.pal(8, "Dark2"))

# barplot
ggplot(d, aes(y = -log10(adj.p.value) , x = reorder(GO.term.name , -adj.p.value) , 
              fill =num_intersection ) ) +
  theme_bw() +
  theme(axis.text = element_text(size = 15) , 
        axis.title = element_text(size = 20), 
        axis.title.y = element_blank() ) + 
  geom_bar(stat = "identity" ) +
  geom_hline(yintercept = -log10(0.05) ,linetype = 'dashed' , size = 1)+
  scale_fill_gradient(low="red", high="blue") +
  coord_flip() 

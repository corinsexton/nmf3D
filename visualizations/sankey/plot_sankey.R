setwd("~/Documents/UNLV/Year4/nmf3D/visualizations/sankey")

library(tidyverse)
library(UpSetR)
library(networkD3)

net <- read_tsv("network.tsv")
network <- as_tibble(table(net))
colnames(network) <- c("source", "target", "value")
network$source <- as.numeric(network$source)
network$target <- as.numeric(network$target)
network$value <- as.numeric(network$value)

# MUST BE ZERO-INDEXED
network$source <- network$source - 1
network$target <- network$target - 1 + 6

# network$target <- network$target + 10

#sm_net <- network[network$n > 200 & network$n < 500,]
# sm_net$target <- sm_net$target + 15

## H3K9ac - active promoters
## H3K27me3 - repressive
## H3K4me3, H3K27me3 - poised promoter
## H3K4me3 - promoter
## H3K27ac - enhancer
## relatively higher H3K4me1 signal compared to H3K4me3 - enhancer


nodes = data.frame("name" = 
                     c(paste0(1:6,""),
                       paste0(1:6,"")),
                   "node" = 0:2,stringsAsFactors = F
)


names(network) = c("source", "target", "value")

# Plot
sn <- sankeyNetwork(Links = data.frame(network), Nodes = nodes,
                    Source = "source", Target = "target",
                    Value = "value", NodeID = "name",fontSize= 14, nodeWidth = 5)

# sn <- htmlwidgets::prependContent(sn, htmltools::tags$h1("HERVH Candidate Loci"))
sn <- htmlwidgets::prependContent(sn, htmltools::HTML('<h1 style="font-family:Arial, sans-serif;font-weight:bold;margin-bottom:0;text-align:center">Label Movement</h1><div style="position:relative;padding: 10px 20px 0px 20px;margin-bottom:50px;">
<p style="font-family:Arial, sans-serif;font-weight:bold;position:absolute;left:0;">H1</p>
<p style="font-family:Arial, sans-serif;font-weight:bold;position:absolute;right:0;">Endoderm</p>
</div>'))

# htmlwidgets::sizingPolicy(padding = 10, browser.fill = TRUE)
# sn$sizingPolicy$viewer$fill <- FALSE
sn$sizingPolicy$browser$fill <- F
sn$sizingPolicy$browser$defaultWidth <- 400

# you save it as an html
saveNetwork(sn, "6lab_nmf.html")

library(webshot)
# you convert it as png
webshot("6lab_nmf.html","6lab_nmf.png", vwidth =400, vheight = 600)


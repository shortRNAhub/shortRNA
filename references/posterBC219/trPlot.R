library(GenomicRanges)

# create object
gr <- new("GRanges",
  seqnames = new("Rle",
    values = structure(1L, .Label = "chr8", class = "factor"),
    lengths = 8L, elementMetadata = NULL, metadata = list()
  ),
  ranges = new("IRanges",
    start = c(
      1L, 1L, 3L, 2L, 51L, 51L,
      55L, 42L
    ), width = c(50L, 36L, 27L, 24L, 28L, 25L, 24L, 19L), NAMES = NULL, elementType = "ANY", elementMetadata = NULL,
    metadata = list()
  ), strand = new("Rle",
    values = structure(3L, .Label = c(
      "+",
      "-", "*"
    ), class = "factor"), lengths = 8L, elementMetadata = NULL,
    metadata = list()
  ), seqinfo = new("Seqinfo",
    seqnames = "chr8",
    seqlengths = NA_integer_, is_circular = NA, genome = NA_character_
  ),
  elementMetadata = new("DataFrame",
    rownames = NULL, nrows = 8L,
    listData = list(abundance = c(
      101, 9982, 59, 310, 109,
      478, 89, 403
    )), elementType = "ANY", elementMetadata = NULL,
    metadata = list()
  ), elementType = "ANY", metadata = list()
)


library(ggplot2)
library(ggbio)
library(cowplot)

theme_set(theme_cowplot())

p1 <- ggplot(gr) + 
  geom_rect(mapping = aes(fill = log(abundance))) + 
  xlab("Position along transcript") + ylab("Unique fragments")

cc <- as.numeric(coverage(gr, weight = gr$abundance)[[1]])
cc <- data.frame(pos = 1:length(cc), Coverage = cc)
p2 <- ggplot(cc, aes(pos, Coverage)) + geom_line(colour = "blue", size = 2)+
  xlab("")


p1 + theme(text = element_text(size=25),
           axis.text.x = element_text(angle=90, hjust=1, size=20), legend.title=element_text(size=20))+
  guides(fill=guide_legend(title = "abundance\n(log)", size = 10))


p2 + theme(text = element_text(size=25),
           axis.text.y = element_text(angle=0, hjust=1, size=20),
           axis.text.x = element_text(angle=0, hjust=1, size=0),
           axis.ticks.x =element_blank(),
           axis.line.x = element_blank()) 

# library(gridExtra)
# 
# grid.arrange(p2, p1@ggplot, nrow = 2, )
# 
# library(ggpubr)
# 
# ggarrange(p2, p1@ggplot + rremove("x.text"), 
#           ncol = 1, nrow = 2, align = "h")
# 
# p1.p <- plotly::ggplotly(p2)

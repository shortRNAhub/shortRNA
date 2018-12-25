library(hexSticker)

imgurl <- "img.png"

cols <- RColorBrewer::brewer.pal(n = 9, name = "Set1")

sticker(imgurl, package="shortRNA", p_size=8, s_x=1, s_y=.8, s_width=0.6, s_height=0.5,
        filename="baseplot.png", h_fill = cols[1], h_color = cols[2])

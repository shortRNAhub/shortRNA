library(hexSticker)

imgurl <- "sRNA.png"

cols <- RColorBrewer::brewer.pal(n = 9, name = "Set1")

sticker(subplot = imgurl, package="shortRNA", p_size=8, p_color = "black", p_family = "Arial", 
        url = "mansuylab/shortRNA",u_color = cols[1], u_family = "Arial",
        s_x=1, s_y=0.8, s_width=0.6, s_height=0.3, h_size = 1.5,
        filename="baseplot.png", h_fill = "white", h_color = cols[2], spotlight = F, u_size = 2)

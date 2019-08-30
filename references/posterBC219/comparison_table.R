file <- readxl::read_excel("./files/Tools_comparison.xlsx", sheet = 2)
DT::datatable(file, rownames = F)

library(kableExtra)
kable(file, format = "html")

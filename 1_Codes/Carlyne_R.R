library(dplyr);library(readr)

f_source <- "E:/DR ASAMOAH/CLIMATE_RISK_WIO/TestData/"
flist <- list.files(f_source, pattern = ".csv")

data1 <- lapply(
  1:length(flist), FUN = function(x){
    as.data.frame(readr::read_csv(paste0(f_source, flist[[x]]),
                           col_types = cols(Study_period_year = col_double())))[,-13]  }) %>% 
  bind_rows()

data1[order(data1$Identi,decreasing = F), ]

#Consider all covers...
data1$totalCover <- apply(data1[,c("urbanCover_10000m","urbanCover_1000m","urbanCover_5000m")],1, FUN = function(x) sum(x,na.rm = TRUE))

# those that duplicated, return 
(data_fnl <- data1 %>% group_by(Identi) %>% arrange(dplyr::desc(totalCover)) %>% filter(row_number()==1) %>% ungroup())
summary(duplicated(data_fnl$Identi))
# Mode   FALSE 
# logical    6063
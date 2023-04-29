#' plot data
#'
#' @description creates bar plot of data profile for each sample, 
#' where x-axis are 96 substitution bases and y-axis are their relative contribution.
#'
#' @param context
#' @param what
#' @param sample_ids
#' @param x basilica object
#'
#' @return plot
#' @export plot_data
#' @examples



plot_data = function(x,sample_ids = NULL,what = "SBS", context = T){
  
  a = x$input$counts
  
  a =  a %>% dplyr::mutate(sbs = rownames(a)) %>% as_tibble() %>%
    reshape2::melt()
  
  
  if(what == 'SBS'){
    a = a %>% dplyr::rename(Var1 = sbs, Var2 = variable) %>%
      mutate(
        substitution = paste0(substr(start = 3, stop = 3, Var2),">",substr(start = 5, stop = 5, Var2)),
        context = paste0(
          substr(start = 1, stop = 1, Var2),
          '_',
          substr(start = 7, stop = 7, Var2)
        )
      )  }
  
  if(what == "DBS"){
    
    a = a %>% dplyr::rename(Var1 = sbs, Var2 = variable) %>%
      mutate(
        substitution = paste0(substr(start = 1, stop = 2, Var2),">NN"),
        context = substr(start = 4, stop = 5, Var2)
      ) 
  }
  
  if(what == "ID"){
    
    a = a %>% dplyr::rename(Var1 = sbs, Var2 = variable) %>%
      mutate( Var2 = as.character(Var2),
              substitution = substr(start = 1, stop = nchar(Var2) - 2, Var2),
              context = substr(start = nchar(Var2), stop = nchar(Var2), Var2)
      ) 
  }
  
  if(what == "CNV"){
    
    a = a %>% dplyr::rename(Var1 = sbs, Var2 = variable) %>% rowwise() %>% 
      mutate( Var2 = as.character(Var2),
              substitution = 
                paste0(str_split(Var2,pattern = ":")[[1]][1],":",str_split(Var2,pattern = ":")[[1]][2]),
              context =  paste0(str_split(Var2,pattern = ":")[[1]][3])
      ) 
  }
  
  
  
  library(CNAqc)
 
  if(!is.null(sample_ids)){ a = a %>% filter(Var1 %in% samples_ids) }
  
  all_plot =  lapply(a$Var1 %>% unique(), function(var){
    
    p =   a  %>% filter(Var1 == var) %>% 
      ggplot() +
      geom_bar(aes(value, x = context, fill = "Indianred"), stat = 'identity') +
      facet_grid(~factor(substitution,levels = a$substitution %>% unique()), scales = 'free') +
      CNAqc:::my_ggplot_theme()  + 
      theme(axis.text.x = element_text(angle = 90)) + 
      guides(fill = 'none')  +
      labs(y = "", title = var)
    
    
    if(!context) {
      p = p + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + labs(x = "")
    }
    
    p
    
  })
  
  all_plot
  
}




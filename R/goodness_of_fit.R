

get_reconstructed_data = function(x){
  
  signatures = NULL
  if("catalogue_signatures" %in% names(x$fit)){ signatures = rbind(signatures,x$fit$catalogue_signatures) }
  if("denovo_signatures" %in% names(x$fit)){ signatures = rbind(signatures,x$fit$denovo_signatures)}
  
  signatures = signatures[colnames(x$fit$exposure),]
  
  as.matrix(x$fit$exposure[rownames(x$input$counts),]*rowSums(x$input$counts)) %*%  as.matrix(signatures) %>%
    as.data.frame()
  
}

get_similarity_scores = function(reconstr_data,real_data){
  
  reconstr_data = reconstr_data[rownames(real_data),]
  cosines = lapply(1:nrow(reconstr_data),function(j){
    
    tibble(id = rownames(real_data)[j], cos =  cosine.vector(reconstr_data[j,],real_data[j,]))
    
  }) %>% bind_rows()
  
  cosines
  
}



plot_gof = function(scores){

if(!"group" %in% (scores %>% colnames())){  scores = scores %>% mutate(group = "1")}
  
  ggplot(scores, aes(x = group,y =  cos,fill = group)) + geom_boxplot() + ylim(0,1) +
   CNAqc:::my_ggplot_theme() + labs(title = "GOF of data counts", y = "GOF")  + ggsci::scale_fill_lancet()
  
}




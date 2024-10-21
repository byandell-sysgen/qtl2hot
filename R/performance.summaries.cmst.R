performance.summaries.cmst <- function(out, model, alpha=0.05, method)
{
  ntests <- nrow(out[[1]])
  if((model=="A") || (model=="B")){
    ct <- counts(out, alpha, method)
    tp <- ct[1,1]
    fp <- ct[1,2]+ct[1,3]+ct[1,4]
    power <- tp/ntests
    type1.err <- fp/ntests
    prec <- tp/(tp+fp)
  } 
  if(model=="C"){
    ct <- counts(out, alpha, method)
    tp <- ct[1,4]
    fp <- ct[1,1]+ct[1,2]+ct[1,3]
    power <- tp/ntests
    type1.err <- fp/ntests
    prec <- tp/(tp+fp) 
  }
  if(model=="D"){
    ct <- counts(out, alpha, method)
    tp <- ct[1,3]
    fp <- ct[1,1]+ct[1,2]+ct[1,4]
    power <- tp/ntests
    type1.err <- fp/ntests
    prec <- tp/(tp+fp)  
  }
  if(model=="E"){
    ct <- counts(out, alpha, method)
    tp <- ct[1,4]
    fp <- ct[1,1]+ct[1,2]+ct[1,3]
    power <- tp/ntests
    type1.err <- fp/ntests
    prec <- tp/(tp+fp)   
  } 
  data.frame(tp, fp, power, type1.err, prec)
}

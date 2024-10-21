performance.summaries.cit <- function(out, model, alpha=0.05)
{
  ntests <- nrow(out[[1]])
  if((model=="A") || (model=="B")){
    ct <- counts(out, alpha, method="cit")
    tp <- ct[1,1]
    fp <- ct[1,2]+ct[1,3]
    power <- tp/ntests
    type1.err <- fp/ntests
    prec <- tp/(tp+fp)
  } 
  if((model=="C") || (model=="D") || (model=="E")){
    ct <- counts(out, alpha, method="cit")
    tp <- ct[1,3]
    fp <- ct[1,1]+ct[1,2]
    power <- tp/ntests
    type1.err <- fp/ntests
    prec <- tp/(tp+fp)   
  }
  data.frame(tp, fp, power, type1.err, prec)
}

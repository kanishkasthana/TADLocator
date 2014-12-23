tailbins=function(stringency){

tailbins=apply(all_rates,2,function(rates){
  mn=mean(rates);
  s=sd(rates);
  new_rows_left[rates<=(mn-stringency*s)]
})
 
return(tailbins[[6]])

}

out=sapply(1:length(sts),function(i)
  {
   
   out=length(intersect(rep1[[i]],rep2[[i]]))/min(length(rep1[[i]]),length(rep2[[i]]));
   
   return(out)
   
  })
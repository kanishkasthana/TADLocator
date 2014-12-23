  
data=read.table('h1.merge.chr17.nor.qq.matrix');

#Removing first three columns and storing it in a new data frame
new_data=data[,c(-1,-2,-3)]

#Converting data to a matrix

interaction_matrix=data.matrix(new_data);
print(dim(interaction_matrix))
orignalDimensions=nrow(interaction_matrix)
#Removing rows with zeros for analysis:

#Getting Rows that have no interaction scores
sums=apply(interaction_matrix,1,sum)
#Removing rows with no interaction and testing again
interaction_matrix=interaction_matrix[!(sums==0),!(sums==0)]
rows_left=(1:orignalDimensions)[!(sums==0)]

print(dim(interaction_matrix))

end=nrow(interaction_matrix)

bias1=numeric(length=end)


for(i in 1:end)
{
  upstream_interaction_count=sum(interaction_matrix[i,1:i])
  downstream_interaction_count=sum(interaction_matrix[i,i:end])
  bias1[i]=upstream_interaction_count-downstream_interaction_count
}

#Computing windowed Bias scores

bias2=numeric(length=end)
#Defining windows for which operation needs to be carried out
windows=c(400,300,200,150,125,100,95,90,80,75,50,25,20,15,10,5)

#This commands returns all the computed rates of change of bias score for all windows mentioned above.
all_rates=sapply(windows, function(window){
  
  print(window)
  for(i in 1:end)
  {
    if(i<=window)
    {
      upstream_interaction_count=sum(interaction_matrix[i,1:i])
      downstream_interaction_count=sum(interaction_matrix[i,i:(i+window)])
      bias2[i]=upstream_interaction_count-downstream_interaction_count
    }
    
    else if(i+window>=end)
    {
      upstream_interaction_count=sum(interaction_matrix[i,(i-window):i])
      downstream_interaction_count=sum(interaction_matrix[i,i:end])
      bias2[i]=upstream_interaction_count-downstream_interaction_count
    }
    else{
      upstream_interaction_count=sum(interaction_matrix[i,(i-window):i])
      downstream_interaction_count=sum(interaction_matrix[i,i:(i+window)])
      bias2[i]=upstream_interaction_count-downstream_interaction_count
    }
  }
  
  rates=diff(bias2);
  return(rates)
  
}
)

#Adding no window results to all_rates to calculate P-value when no window is used:
all_rates=cbind(all_rates,diff(bias1))

new_rows_left=rows_left[2:length(rows_left)]

sts=seq(2,2.5,by=0.05)

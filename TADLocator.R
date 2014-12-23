#Script for Analysing and locating TADs by Kanishk Asthana
#Reading in Interaction File for Chromosome 17
require('MESS')
data=read.table('h1.merge.chr17.nor.qq.matrix');

#Removing first three columns and storing it in a new data frame
new_data=data[,c(-1,-2,-3)]

#Converting data to a matrix

interaction_matrix=data.matrix(new_data);

pdf("Rep1.pdf")

#interaction_heatmap <- heatmap(log2(interaction_matrix), Rowv=NA, Colv=NA, col = rev(heat.colors(256)),scale="none",revC=TRUE, main="HeatMap of Interactions")


print(dim(interaction_matrix))
orignalDimensions=nrow(interaction_matrix)

#Removing rows with zeros for analysis:
#Getting Rows that have no interaction scores
sums=apply(interaction_matrix,1,sum)
#Removing rows with no interaction and testing again
interaction_matrix=interaction_matrix[!(sums==0),!(sums==0)]
rows_left=(1:orignalDimensions)[!(sums==0)]

print(dim(interaction_matrix))

#scaled_heatmap <- heatmap(log2(interaction_matrix), Rowv=NA, Colv=NA, col = rev(heat.colors(256)),scale="none",revC=TRUE, main="Heatmap after removing missing data")

print(dim(interaction_matrix))

end=nrow(interaction_matrix)

bias1=numeric(length=end)


for(i in 1:end)
{
   upstream_interaction_count=sum(interaction_matrix[i,1:i])
   downstream_interaction_count=sum(interaction_matrix[i,i:end])
   bias1[i]=upstream_interaction_count-downstream_interaction_count
}


hist(diff(bias1),500,main='Histogram of Rates with no windowing',xlab='Rates')

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
 
#Plotting Histogram for each window

rates=diff(bias2);
test=shapiro.test(rates)
print(test$statistic)

#hist(rates,250,main=paste('Histogram of Rates with window of ',window),xlab='Rates')
return(rates)

}
)

#Adding no window results to all_rates to calculate P-value when no window is used:
all_rates=cbind(all_rates,diff(bias1))

#Getting all domain boundaries:
boundaries=read.table('finaldomainclls.chr17');
#Storing all boundaries as a matrix:
boundaries_matrix=as.matrix(boundaries[-1])

#Getting location of bins at boundaries
boundary_bins=boundaries_matrix/40000

#Getting unique position of Bins in this distribution
unique_bins=sort(unique(boundary_bins[1:(2*nrow(boundary_bins))]))


#Defining function to check if bin is present in already discovered TAD boundaries.
#Here a bin is said to be present at the boundary if the bin is present at +1 or -1 bins of a discovered TAD boundary.


#Correcting for the fact that bins begin at Zero in the Data File
unique_bins=unique_bins+1

isinbins= function(bin){
	if( ( sum(unique_bins==bin) + sum(unique_bins==(bin+1)) + sum(unique_bins==(bin-1)) ) > 0 ){
	return(1)
	}
	else
	{
	return(0)
	}
}

#Defining function that gives us the fraction of an array of bins present at the boundaries

fractioninboundaries= function(bin_vector){

nelems=length(bin_vector)

inbins=sapply(bin_vector,isinbins)

fraction=sum(inbins)/nelems

return(fraction)

}


#Getting number of samples which are less than -100 in rates
bootstrapsizes=apply(all_rates,2,function(rates){ mn=mean(rates);s=sd(rates);sum(rates<=(mn-3.5*s))})

print(bootstrapsizes)

#Making bootstrap confidence intervals for the given sample sizes:
Bsims=1000

#Calculating Bootstrap confidence intervals
confintervals=sapply(bootstrapsizes,function(size){

resamples= matrix(sample(rows_left,Bsims*size,replace=TRUE),Bsims,size)
print(dim(resamples))

resamplefractions=apply(resamples,1,fractioninboundaries)
print(length(resamplefractions))

return(resamplefractions)

}
)

apply(confintervals,2,function(intervals){ plot(density(intervals))})

#Calculating fraction TADs in the tails of the distribution:

#Correcting for the fact that taking a diff reduces the length of the vector by 1
new_rows_left=rows_left[2:length(rows_left)]

tailbins=apply(all_rates,2,function(rates){
  mn=mean(rates);
  s=sd(rates); 
  new_rows_left[rates<=(mn-3.5*s)]
})

tailsinboundaries=sapply(tailbins,fractioninboundaries)

print('Percentage located at already known tail boundaries:')
print(tailsinboundaries)

#Now calculating P-values: 

pvalues=sapply(1:length(tailsinboundaries),function(i){
#Value calculated for which p-value need to be calculated from bootstrap confidence intervals
val=tailsinboundaries[i]
return(sum(confintervals[,i]>=val)/Bsims)
})


print('The P values for all the windows defined above are:')
print(pvalues)



#Parsing out the two distributions to generate ROC curve

vector=logical(length=length(rows_left))

for(i in 1:length(unique_bins))
{
 vector=vector+(rows_left==unique_bins[i])
}

#Vector now contains an array of elements with the bins present in the cleaned data

#Logical Vector of values
logical_vector=vector>0

#Creating ROC Curves

#Modifying Logical Array such that we get one element after and one element before each boundary
i=0

repeat{ 
i=i+1;

 if(logical_vector[i]==TRUE){
  if(i!=1){ logical_vector[(i-1)]=TRUE }
  if(i!=length(logical_vector)){  logical_vector[i+1]=TRUE }
  i=i+2
 }
 if(i==length(logical_vector)){break}
}

inverse_logical=!logical_vector
#Shortening for diff
logical_vector=logical_vector[2:length(logical_vector)]

inverse_logical=inverse_logical[2:length(inverse_logical)]

#Creating range of values for which ROC curve is to be made.
#This was a lot of fun to do, it involved a lot of problem solving, that I enjoyed.
#Well the point it to have fun.


divisions=seq(-300,300,by=1)

rocs=apply(all_rates,2,function(rates){

#Getting Distribution of Rates for Boundaries 
truedist=rates[logical_vector]

#Getting Distribution of rates Non-boundaries
dist=rates[inverse_logical]

#Getting Total Number of Elements in the Distribtuion of Boundaries
allnums=sum(logical_vector)

#Getting Total Number of Elements in Non-boundaries
distnums=sum(inverse_logical)

#Getting a vector for sensitivity
sensitivity=sapply(divisions,function(div){

return(sum(truedist<=div)/allnums)
 
})
#Getting a vector for specificity
specificity=sapply(divisions,function(div){
return(sum(dist>div)/distnums)

})
#Plotting ROC curve for given Rates in a given window: No window case is the last one
plot(1-sensitivity,specificity)

print(auc(1-sensitivity,specificity,type='spline'))

#Returning sensititvity and specificity values
return(cbind(sensitivity,specificity))
})


dev.off()


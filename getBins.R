#Function created by Kanishk Asthana
# Defining Function getBins to get the Bin numbers for Domain boundaries in a MATRIX file for a specific stringency level
#Here stringency is defined as the number of standard deviation from the mean of the distribution of rates of the windowed Bias Score
#So for example a stringency of 1 will give you all the bins that have a rate lower than -1*(standard deviation) of the histogram
#of rates described in the report


getBins <- function(filename,stringency=2.5,window=100){
  #Reading file
  data=read.table(filename);
  
  #Removing first three columns and storing it in a new data frame
  new_data=data[,c(-1,-2,-3)]
  
  #Converting data to a matrix
  #All operations done from now onwards assume that we are dealing with a square matrix
  interaction_matrix=data.matrix(new_data);
  
  orignalDimensions=nrow(interaction_matrix)
  
  #Removing rows with zeros for analysis:
  
  #Getting Rows that have no interaction scores
  sums=apply(interaction_matrix,1,sum)
  #Removing rows and columns with no interaction and testing again
  interaction_matrix=interaction_matrix[!(sums==0),!(sums==0)]
  #Storing Row numbers that were not deleted
  rows_left=(1:orignalDimensions)[!(sums==0)]
    
  #Storing number of rows and columns of the matrix assuming it is a square matrix.
  end=nrow(interaction_matrix)
  
  #Computing windowed Bias scores
  bias2=numeric(length=end)
  
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
  
  
  #Computing Rates of Change of Scores for given window
  rates=diff(bias2);
  
  #Correcting for the fact that taking a diff reduces the length of the vector by 1
  new_rows_left=rows_left[2:length(rows_left)]
  
  #Getting the bins for given stringency
  mn=mean(rates);
  s=sd(rates); 
  new_rows_left[rates<=(mn-stringency*s)]
  
}
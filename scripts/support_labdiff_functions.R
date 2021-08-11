average_labelling_proportion <- function(df, targets)
{
  avg_lab_df <- sapply(unique(targets[,2]), function(x, df, targets){
    samples <- targets[targets[,2] == x,1]
    sub_df <- df[,samples]
    avg <- rowMeans(sub_df)
    return(avg)
  },df = df, targets = targets)
  
  return(avg_lab_df)
}

average_labelling_differences <- function(avg_lab_df, comparisons)
{
  avg_diff_df <- lapply(comparisons, function(x, avg_lab_df){
    avg_diff <- avg_lab_df[,1, drop = F]
    avg_diff[,1] <- 0 
    for(index in x)
    {
      # print(avg_diff)
      # print(avg_lab_df[,index, drop = F])
      # print(index)
      if(index > 0)
      {
        avg_diff <- avg_diff + avg_lab_df[,index, drop = F]
      } else
      {
        avg_diff <- avg_diff - avg_lab_df[,abs(index), drop = F]
      }
    }
    # View(avg_diff)
    return(avg_diff)
  }, avg_lab_df = avg_lab_df)
  return(avg_diff_df)
}
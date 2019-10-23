F_JointPvalueFromObserved <- function(testTab){
   # testTab must contain
   # a column $tau giving the probabilty to jointly H1 for the pair
   # a column $stat giving the joist statistic
   # will return the same table + a column $pval
   # testTab = jointTests
   
   nbTests = nrow(testTab); testTab$rowNum = (1:nbTests)
   testTab = testTab[order(testTab$stat, decreasing=TRUE), ]
   testTab$pval = (.5 + (1:nbTests) - cumsum(testTab$tau)) / (1 + nbTests - sum(testTab$tau))
   testTab = testTab[order(testTab$rowNum), ]; testTab$rowNum = c()

   return(testTab)
}


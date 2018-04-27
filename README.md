1. Description: 
   Detecting differentially expressed genes.

2. Usage:
(1). Install packages: "edgeR" and "fitdistrplus" in R (>=3.2.0).
(2). Run "q.value" function.
(3). Call DESyn(data, group).

3. Arguments of DESyn:
   (1). data: a data frame containing gene IDs and normalized counts. The first column gives gene ID's, and the other columns provide normalized counts.
   (2). group: vector of group indicators.
   
   Example:
   (1). data  
   ID     Sample1     Sample2     Sample3     Sample4     Sample5     Sample6
   g1     17          21          31          15          50          12
   g2     89          10          79          109         70          28
   (2). group=c(0,0,0,1,1,1) #group indicator equals to 0 if the sample is from the control group and 1 if the sample is from the case group.

   Note: We recommend to have at least 3 samples from each group.

4. Value: 
   (1). geneID: a vector containing geneID in the dataset.
   (2). pvalue: p-value of genes.
   (3). qvalue: q-value of genes.


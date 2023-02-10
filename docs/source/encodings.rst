Encodings
---------

To use Bayesian optimization, LC/MS data must be encoded. Here, we encode the LC gradient as 
the input to the model, and the MS data as the output. 

The LC gradient is encoded as a numeric vector, with each element in the vector representing
the mobile phase ratio at a given time point. 

The MS data were processed and encoded as a singular value between 0-1, where 0 represents
no separation and 1 represents complete separation.

Read more details about the encoding process in the paper:



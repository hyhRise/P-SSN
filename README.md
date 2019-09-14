# P-SSN
 Construct the single-sample network based on partial correlation(P-SSN)

## construct reference network
The "construct_reference_network.py" can be used to construct the background network,and the parameter is shown below:

“construct_reference_network.py” requires "numpy" extensive package.

usage: python construct_reference_network.py -process=process_value -proportion=proportion_value -PCC=threshold_PCC -ref=reference_sample_file  -sample=sample_data_file -out=results_output_fold
Options and arguments:
-process: namesumber of work processes used
-proportion: the genes whose expression value was 0 in more than proportion samples/cells
-PCC : set the threshold of PCC [0..1], if the -pvalue set 0, all edges will be outputted to the background_network
-ref : the expression profile of reference samples
-sample : the expression profile for the sample to be constructed the PSSN
-out : the file to store the background_network

There is a simple example below: 
The file "reference_samples.txt" is the expression profile of reference samples.
The file "14_samples_data.txt" is a 14 samples expression profile.


For example,to construct the P-SSN for the 14 samples, and the results will be put in "reference_network.txt" file:

python construct_reference_network.py -process=2 -proportion=0.5 -PCC=0.7 -ref=reference_samples.txt  -sample=14_samples_data.txt -out=reference_network.txt


----------------------------------------------------------------------------------------------------------------------------------------

## construct single network
The "construct_single_network.py" can be used to construct the single-sample network based on partial correlation(P-SSN), and the parameter is shown below:

The "construct_single_network.py" requires "numpy" and "scipy" extensive package.

usage: python construct_single_network.py -proportion=proportion_value -pvalue=threshold_p_value -background=background_network_file -ref=reference_sample_file  -sample=sample_data_file -out=results_output_fold
Options and arguments:
-proportion: the genes whose expression value was 0 in more than proportion samples/cells
-pvalue : set the threshold of p-value [0..1], if the -pvalue set 1, all edges will be outputted to the P-SSN
-background : background network to calculate the deltaPTCC of edges based on the network
-ref : the expression profile of reference samples
-sample : the expression profile for the sample to be constructed the P-SSN
-out : the directory to store the P-SSN

There is a simple example below: 
The file "reference_samples.txt" is the expression profile of reference samples.
The file "background.txt" is the background network for calculating P-SSN.
The file "14_samples_data.txt" is a 14 samples expression profile, the P-SSNs will be constructed for the 14 samples profile.


For example,to construct the P-SSN for the 14 samples, and the results will be put in "target" directory:

python construct_single_network.py -proportion=0.5 -process=2 -ref=reference_samples.txt -background=reference_network.txt -sample=14_samples_data.txt -pvalue=0.01 -out=target


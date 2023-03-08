# Evaluation-of-FMD-control-measures-in-Thailand
The code for FMD outbreak simulation in Thailand.  
The data that support the findings of this study are available on request from the corresponding author.  
Thanicha Chanchaidechachai Thanicha.c@outlook.co.th

Supplementary S1 - S4
Supplementary S1 Outbreak data
	The outbreak data were collected from 2 areas of Thailand with the FMD outbreak between 2016 and 2017. The first area is Lamphaya Klang subdistrict. The size of the study area in Lamphaya Klang subdistrict is 12.5 × 8.4 km2 covering 502 dairy farms. The FMD outbreak happened from 15 September 2016 to 8 August 2017. The infected farms included 273 dairy farms. The second area is Bo Phloi district. the size of the study area in Bo Phloi district is 30.8 × 25.5 km2 covering 346 beef cattle farms, 104 goat farms and 51 pig farms. The FMD outbreak happened from 13 October 2016 to 15 December 2016. The detail data The infected farms included 15 beef farms. The detail on the questionnaires and infection date estimation can be found in https://doi.org/10.1016/j.prevetmed.2021.105468.  The raw data of farm type, farm size and outbreak date is available for download in the file S1.1LP_outbreak for Lamphaya Klang subdistrict and S1.2BP_outbreak for Bo Phloi district. Due to confidential data, we cannot provide the farm location. We only provided the matrix of between-farm distance (metres unit) in file for S1.3LP_distance for Lamphaya Klang subdistrict and S1.4BP_distance for Bo Phloi district.
S2 Trade network
	The raw data of trade network can be found in file S2.1Trade_network. Due to confidential data, we used ID to identify the traders. Overall, 64 farms traded and could recall the traders' names in Lamphaya Klang subdistrict. There are 48 distinct traders in the area. For Bo Phloi district, 32 farms traded and could recall the traders' names. There are 16 distinct traders.  
Regarding the trade network simulation, we simulated the trade networks under the assumption that 50% of the farms traded with at least one of the traders from the interview. We assumed that farms could only have edges with traders, and the number of edges of each farm was randomized from the degree distribution of the trade network from interview data. The edges between farms and traders were randomly generated using the preferential attachment. The code for trade network simulation can be found in https://github.com/AnnThanicha/Evaluation-of-FMD-control-measures-in-Thailand/blob/main/S2%20trade%20network%20simulation.R
S3 Transmission parameters estimation
	The transmission parameter were estimated as the methods explained in our previous paper https://doi.org/10.1016/j.prevetmed.2021.105468.  To estimate the transmission parameters, the file for outbreak data, trade network and distance between farms are required. The code for parameter estimation and sample data can be found in https://github.com/AnnThanicha/Evaluation-of-FMD-control-measures-in-Thailand/blob/main/S3%20transmission%20parameters%20estimation.R
S4 R Model simulation
	Before running the R code for simulation, the farm input data, trade matrix and between farm distance matrix need to be loaded.
The R code for baseline simulation can be found in https://github.com/AnnThanicha/Evaluation-of-FMD-control-measures-in-Thailand/blob/main/S4%20Baseline%20simulation.R
The R code for culling animals in infected farms can be found in https://github.com/AnnThanicha/Evaluation-of-FMD-control-measures-in-Thailand/blob/main/S4%20Culling%20simulation.R
The R code for ring vaccination can be found in 
https://github.com/AnnThanicha/Evaluation-of-FMD-control-measures-in-Thailand/blob/main/S4%20Ring%20vaccination%20simulation.R
The R code for animal movement restrictions be found in https://github.com/AnnThanicha/Evaluation-of-FMD-control-measures-in-Thailand/blob/main/S4%20Animal%20movement%20restrictions%20simulation.R
The R code for isolation of infected farms be found in
https://github.com/AnnThanicha/Evaluation-of-FMD-control-measures-in-Thailand/blob/main/S4%20Isolation%20simulation.R

 




# PCMLimbLength

Scripts and data to run analyses on hindlimb length allometry.

Scripts PCMFit_MixedModels.R and MixedModelSensitivity.R were run on an R server in parallel running R version 4.3.1.
All other scripts were run on a PC (even if in parallel) running R version 4.3.2.
Data files include: 
1. final_mm.RDS (individual species measurements)
2. ext2_td.RDS (a treeplyr tree data object which combines a phylogeny from Title et al. 2024, data from Meiri 2018, and species means from our dataset)
3. All other .RDS files were generated through the use of the below scripts. 

Scripts should be run in this order:
1. MixedModels/PCMFit_MixedModels.R
2. MixedModels/MixedModelsSensitivity.R
4. MixedModels/fits&plots.R
5. HabitatModels/substrateSIMMAP.R
6. HabitatModels/fits&plots.R

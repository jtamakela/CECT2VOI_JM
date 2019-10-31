# CECT2VOI_JM

Matlab code for analyzing registered triple-contrast agent micro-CT images

The order of execution:
1. registration_check_JM.m 
   - Plots mid slices of the .nii files. Can be used to check the registrations
2. CECT2VOI.m 
   - Takes the predetermined 6 locations, and with the help of the user, creates rotated volumes of interest
3. CECTVOI_Rotate.m 
   - Further rotates the previously created volumes of interest to make sure that the sample is straight. Uses CECTVOI_Alignment.m and find_demarcation.m 
4. VOI_check_JM 
   - Use this to check the VOIs. Displays 'em
5. CECTVOI_Analysis.m 
   - Is used to analyse the previously created VOIs and create depth-dependent attenuation profiles
6. Partition_Analysis.m
   - Takes the depth-dependent profiles and calculates partition profiles

Files with 'check' are not necessary in the analysis. 
# Fellows_2023
All scripts used in Fellows et al 2023 to analyse photobleaching steps
This contains 5 scripts all used for spot intensity measurements. To run the analysis you will require a .tif stack of your image and a separate .csv file with the x/y coordinates of the spots you wish to measure. These should be in the same folder and have the simialer name. E.g video_1.tif, Spots#video_1.csv

All scripts require Matlab to run.

Bleach1_TifSplit - creates separate tif stacks of individual spots.
Bleach3_Measureintensity_Gaussian - measures the intensity of spots. Matlab script reads X and Y values from column 5/6 in the excel file.


Must be present in the same folder as the data:
Bleach3_StepFindIntensity - Calculates changes in intensity and indentifies steps.
Alex4GUI - GUI to visualise spots and step intesnity changes that occur over the tif stack
Alex5Harvest - Collects data from individual spots and collates into one file.



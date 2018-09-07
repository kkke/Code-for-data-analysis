# Code-for-Imaging data-analysis
Codes here are used for analyzing imaging data recorded with Mscan and Intan.

Specifically, all Imaging files and video files (*.MDF) are recorded with MScan at 30 Hz, and the behavioral event such as licking, taste delivery and scanning signals are recorded with Intan (*.rhd).

To analyze the data,

First, you need to export Imaging file (*.MDF) to TIFF files with MView, rename the first 9 episodes of imaging files as *01.tif,*02.tif,.. and *09.tif instead of *1.tif, *2.tif,...and *9.tif;

Second, launch FIJI, go the macro menu and install the macro called batchAverag.ijm. Launch the batchAverag from the macro list,and windows will pop out to guide you to automatically temporally downsample the tiff files by group averaging (5 times donwsamplling). Save the output downsampled tiff files.

Third, run the suite2P package to perform autmatically segmentation: segmentation_suit2p.m

Forth, run new_main.m to manually correct the segmentation.

Fifth, run imaging_analysis_GC_v4.m to analyze data for GC project.



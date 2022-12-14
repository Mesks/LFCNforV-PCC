# Lightweight fully connected network-based fast CU size decision for video-based point cloud compression
This is the official repository of source codes and deployment methods for the paper "Lightweight fully connected network-based fast CU size decision for video-based point cloud compression". In order to reduce its meaning and express it uniformly, the following takes "LFCNforV-PCC" as the root directory, for example, the location of "mpeg-pcc-tmc2-netTest" is "/mpeg-pcc-tmc2-netTest"

<b>If you have contacted TMC2 and HM, you can skip subsequent lengthy instructions and directly use the source under “/mpeg-pcc-tmc2-netTest/dependencies/HM-16.20+SCM-8.8/source” to change the source of original dependencies of TMC2 and check the methods described in our paper.  If you are not familiar with the package structure of TMC2, it is strongly recommended that you configure it as described below.</b>

<b>Or if you just want to run the program to verify the experimental data in the paper, firstly, unzip "pointCloudCompress_testSequence.zip" to the same location as the two empty folders in order to ensure that the batch files run properly, and then you can run "/__batchProcessing/_BP_1.bat" and "/__batchProcessing/_BP_2.bat" directly after configuring the input point cloud file correctly. In order to properly configure the input point cloud, you need to place the five point cloud sequences provided by MPEG under "/pointCloudCompress_testSequence", and under each sequence's Ply subfolder, place a subfolder named sequencename+_n, For example, the full path of the normal file of the first frame of the soldier sequence is "/pointCloudCompress_testSequence/soldier/soldier_n/soldier_vox10_0536_n.ply"</b>

## <b>Resource Link
The program versions used in the experiment are as follows (You can get their official versions through the link after quotation marks): 

1. TMC2-v18.0: https://github.com/MPEGGroup/mpeg-pcc-tmc2/tree/release-v18.0
2. HDRTools-v0.18: https://gitlab.com/standards/HDRTools/-/tree/0.18-dev
3. MPEG test sequence: https://mpeg-pcc.org/index.php/pcc-content-database/

## <b>Content Introduction
<b>In order to reduce the size of GitHub uploaded files, only some key files are uploaded. And first, unzip "pointCloudCompress_testSequence.zip" to the same location as the two empty folders in order to ensure that the batch files run properly.</b> A brief introduction to the content provided is listed below:  

- __batchProcessing: Store batch files of .sh and .bat, which you can run after configuring the input point cloud to check the experimental data in the paper or configure your own caller using this as a reference.

- __output: Store the intermediate files generated by TMC2, where you can check the occupancy map, geometry map and attribute map. It is necessary to keep the intermediate files. If you want to save space by not saving the intermediate files, you need to set "--keepIntermediateFiles=0" in the configuration file and close the part of the program that gets the placeholder image (the latter is added to the algorithm by us middle).

- __statisticData: The console output is stored, and we provide the statistical tool "__TMC2_statisticDataTool.exe" and the configuration file "__setting.init" based on the console output. Please modify the input file name and output file name in the same directory in the configuration file before use, and normally configure the number of QPs and frames contained in a single file.

- external: Store the external dependencies of TMC2 on HDRTools, you can also download it yourself through the link above.

- LFCN_training: Store the content related to neural network training, including datasets, trained models and history and python files used for training. We provide two versions: py and ipy.

- mpeg-pcc-tmc2-featuresExtracting, mpeg-pcc-tmc2-netTest and mpeg-pcc-tmc2-Xiong: They are the TMC2 program for data extraction, the TMC2 program for neural network testing, and the control experiment Xiong et al.'s algorithm reproduced on TMC2_v18.0. You can directly replace the source file with the source file that HM depends on under dependencies and generate and view it, we did not upload the complete TMC2 program for the sake of reducing the upload volume. In the source file we provide, you can locate the section we have changed by searching for "MesksCode".

- pointCloudCompress_testSequence: If you need to run our batch file directly, you must store the five decompressed sequence files and the corresponding normal files in this folder, otherwise you may not be able to verify the experimental data in the paper.

## <b>Input File
You can download the official test sequences provided by MPEG in the 3 resources mentioned above, but for some reasons, they no longer provide the complete sequence of loot, queen, redandblack, solder, and longdress described in the paper, but basketball_player is still provided. If needed you can still download the sequences mentioned in the paper here: http://plenodb.jpeg.org/pc/8ilabs/

Please decompress the obtained test sequence to the "/pointCloudCompress_testSequence" path. For example, for the input point cloud file of the first frame of the soldier sequence, it should be found by "/pointCloudCompress_testSequence/soldier/Ply/soldier_vox10_0536.ply", and the corresponding The normal file should be able to be found via "pointCloudCompress_testSequence/soldier/soldier_n/soldier_vox10_0536_n.ply". If you are using the normal calculated by MeshLab based on the input file, although it can be run, the compressed PSNR may be different from the PSNR recorded in the official test file, and may not be consistent with our experimental results.

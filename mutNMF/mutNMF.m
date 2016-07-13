cd /data/home/zuoxy/data/NPC/somatic/20160519_863closing/5-mutation_signature/
clear all;
addpath('/data/home/zuoxy/data/test/code/source/');
addpath('/data/home/zuoxy/data/test/code/plotting/');
clc;
allInputFile='1-input/NPC_input.mat';

%% import data
cancerType='NPC';
originalGenomes=csvread('1-input/originalGenomes.csv',1);
subtypes=textread('1-input/subtypes.txt','%s');
types=textread('1-input/types.txt','%s');
sampleNames=textread('1-input/sampleNames.txt','%s');
save(allInputFile)

%% Open matlabpool
if ( matlabpool('size') == 0 )
    matlabpool open; % opens the default matlabpool, if it is not already opened
end

%% Define parameters
iterationsPerCore = 100;
minNumberOfSignature = 1;
maxNumberOfSignature = 10;
stability = zeros(maxNumberOfSignature, 1);
reconstructionError = zeros(maxNumberOfSignature, 1);
inputFile = allInputFile;
allOutputFile = '2-output/res_2_NPC.mat';

%% Sequentially deciphering signatures between minNumberOfSignature and maxNumberOfSignature
for totalSignatures = minNumberOfSignature : maxNumberOfSignature
    
    % Decipher the signatures of mutational processes from catalogues of mutations
    [input allProcesses allExposures idx processes exposures processStab processStabAvg] = ...
        decipherMutationalProcesses(iterationsPerCore, totalSignatures, inputFile, ...
            ['2-output/res_2_NPC_finding_' num2str(totalSignatures) '_signatures.mat'] );
    
    % Record the stability and average Frobenius reconstruction error
    stability(totalSignatures-minNumberOfSignature+1) = mean(processStabAvg);
    reconstructionError(totalSignatures-minNumberOfSignature+1) = norm(input.originalGenomes - processes*exposures, 'fro');
    
end

%% Plotting the stability and average Frobenius reconstruction error
try %% Some versions of MATLAB plotyy has a bug under linux with -nodisplay -nosplash -nodesktop options
  plotSignatureStabilityAndReconstruction(minNumberOfSignature:maxNumberOfSignature, stability, reconstructionError, input);
catch ME
  %% Do not do anything - just ignore the plot in order to save the final output daya
end
%% Saving the data
save(allOutputFile);

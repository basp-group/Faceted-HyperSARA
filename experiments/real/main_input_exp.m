clear ; clc; close all; delete(gcp('nocreate'))
%% change path if needed
hardware = 0 ; %1 if cirrus and 0 if local machine
main_dir='/Users/ad33/CodesScience/data_interface/';
cubedata_filename=[main_dir,filesep,'copyfhs',filesep, 'data',filesep];
calib_dir=[cubedata_filename,'pre_processing_die/'];

% job_dir = [main_dir,'Faceted-Hyper-SARA-EXP2/experiments/real'];
% cd(job_dir)
%% general params
% data 
imagecubeName = 'CygA';
datasetNames={'CYGA-ConfigA','CYGA-ConfigC'};
% data organisation
% option 1: provide 'effChans2Image' a cell array containing the ids of the  channels to be concatenated for each effective channel. 
% option 2: provide all ids of channels 'nChannelsPerImage' & num of channel per effective channel 'nChannelsPerImage' channel
idChannels2Image  = [1:16]; % ids of the channels to be imaged 
nChannelsPerImage   = 16; % number of consecutive channels to be concatenated into each effective channel
nEffChans2Image=floor(numel(idChannels2Image)/nChannelsPerImage);
effChans2Image=cell(nEffChans2Image,1);
for iEff =1:nEffChans2Image
    if iEff<nEffChans2Image, effChans2Image{iEff}=idChannels2Image((iEff-1)*nChannelsPerImage +1:iEff*nChannelsPerImage);
    else,effChans2Image{iEff}= idChannels2Image((iEff-1)*nChannelsPerImage +1:end);
    end
end

dataFilename = @(idSet,ch) strcat(cubedata_filename,filesep,datasetNames{idSet},filesep,'data_ch_',num2str(ch),'.mat');

% running one subcube at a time
subcubeInd=1; % 

% measurement op. params
param_global.measop_flag_wproj =0;
param_global.measop_flag_dataReduction =1;  
%
% image details, dims &  cellsize
param_global.im_Nx = 2560;
param_global.im_Ny = 1536;
param_global.im_pixelSize = 0.06;% pixelsize in asec
%
%faceting params
param_global.facet_Qx =1; % dimFacet1
param_global.facet_Qy =1; % dimFacet2
param_global.facet_Qc =1; % nSubCubes
param_global.facet_window_type='triangular';
param_global.facet_overlap_fraction =[0.1,0.1];
%
% reg params
param_global.reg_gam=0.33; % l21 reg param
param_global.reg_gam_bar=0.33; % nuclear norm reg param
param_global.reg_nReweights=5; % number of re-weights
param_global.reg_flag_reweighting=1; % flag re-weighitng
param_global.reg_flag_homotopy = 0  ;
%
% algo & parallelisation params
param_global.algo_version='sara';
param_global.algo_flag_computeOperatorNorm = 1  ;
param_global.algo_flag_solveMinimization = 1  ;
%
% filenames and input
param_global.exp_type = 'exp2_interface';
param_global.main_dir =[main_dir,'/copyfhs/'];
param_global.preproc_filename_die = @(firstch,lastch) strcat(calib_dir,'/dies/chs',num2str(firstch),'-',num2str(lastch),'_dies.mat');
param_global.preproc_filename_l2bounds = @(firstch,lastch) strcat(calib_dir,'/l2bounds/chs',num2str(firstch),'-',num2str(lastch),'_l2bounds.mat');
param_global.preproc_filename_model = @(firstch,lastch) strcat(calib_dir,'/model_images/chs',num2str(firstch),'-',num2str(lastch),'_model_image.fits');
% hardware
param_global.hardware='local';% 'cirrus' or 'local', add your own cluster


%% run main job
main_real_data_exp(imagecubeName,datasetNames,dataFilename,subcubeInd,effChans2Image,param_global) ;    

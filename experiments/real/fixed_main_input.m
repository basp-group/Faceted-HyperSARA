%% change path 
cirrus = 1;
calib_dir = '/lustre/home/shared/sc004/FACETED_HYPERSARA_EXPERIMENTS/calib/'
cube_filename = '/lustre/home/shared/sc004/FACETED_HYPERSARA_EXPERIMENTS/data/cyga_data_subcube_'

project_dir = pwd;

%% general params
% image details, dims &  cellsize
param_global.Nx = 2560;
param_global.Ny = 1536;
param_global.pixelSize = 0.06 ;%asec
%faceting params
param_global.Qx =5;
param_global.Qy =3;
param_global.Qc =1;
param_global.window_type="triangular";
param_global.overlap_fraction =[0.5,0.5];
%
param_global.exp_type = 'realexp';
param_global.main_dir ='../../';
param_global.G_filename = @(subcube,ch) strcat(calib_dir,'/Gw/subcube',num2str(subcube),'_ch',num2str(ch),'_Gw.mat');
param_global.l2bounds_filename = @(subcube,ch) strcat(calib_dir,'/l2bounds/subcube',num2str(subcube),'_ch',num2str(ch),'_l2bounds.mat');
param_global.model_filename = @(subcube,ch) strcat(calib_dir,'/model_image/subcube',num2str(subcube),'_ch',num2str(ch),'_model_image.fits')


%% flags 
input_flags.computeOperatorNorm = 1  ;
input_flags.solveMinimization = 1  ;
input_flags.homotopy = 0  ;
input_flags.dr =0;
input_flags.wprojection =0


% cd ([project_dir,filesep,'experiments',filesep,'real']);
imagename = 'CygA'  ;

algoversion = 'fhs';

ncoredata = 30  ; %30 for cyga

channels2image = [1:17 19 21:32]; %ids of the channels to be imaged from subcube=subcube_ind;
%if 'sara', make sure to select the id of the channel to be imaged.
ncoredata =param_global.Qx *param_global.Qy+ numel(channels2image)+1; %30 for cyga




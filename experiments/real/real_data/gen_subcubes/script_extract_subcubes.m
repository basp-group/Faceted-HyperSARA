clear; clc; close all;
SpeedOfLight = 299792458;

obsFileName = {'5','7'};
configFileName = {'A','C'};
pixelSize = 0.06;%asec
Qc=16;
direname=['/lustre/home/shared/sc004/FACETED_HYPERSARA_EXPERIMENTS/data',filesep];
file2save = [direname,'cyga_data_subcube_'];
file2load = cell(length(configFileName),length(obsFileName));
for conf = 1 : length(configFileName) %visibility file
    for m = 1 : length(obsFileName)  %visibility file
        file2load{conf,m}=([direname,'CYG-' configFileName{conf} '-' obsFileName{m} '-FULL.mat']);
    end
end
% % util_get_real_vis_data_ad(file2load,file2save,pixelSize,Qc);


for conf = 1 : length(configFileName) %visibility file
    ch =1;
    for m = 1 : length(obsFileName)  %visibility file
        load(file2load{conf,m});
        %% Split data
        % number of overall channels
        c = size(vis,2);
        % number of spectral windows
        s = length(unique(data_id));
        sources = unique(field_id);
        cyg = sources(end);
        
        vis = vis(field_id==cyg,:);
        flaging = flaging(field_id==cyg,:);
        flag_row = flag_row(:,field_id==cyg)';
        uvw = uvw(field_id==cyg,:);
        weights_ch = weights_ch(field_id==cyg,:);
        times = times(:,field_id==cyg)';
        data_id = data_id(:,field_id==cyg)';
        
        for i = 1 : s
            Data = vis(data_id==i-1,:);
            Flagging = flaging(data_id==i-1,:);
            Flag_row = flag_row(data_id==i-1,:);
            UVW = uvw(data_id==i-1,:);
            Weights_ch =(weights_ch(data_id==i-1,:,:));
            Times = times(data_id==i-1,:);
            
            for j = 1:c %ch_ind
                wavelength = SpeedOfLight/Freqs(i,j);
                temp = Flagging(:,j) | Flag_row;
                I = double(Data(:,j));
                if strcmp(obsFileName{m},'7') && strcmp(configFileName{conf},'A') && i ==15
		    ind =(~temp);
		    time_{ch}{conf} = double(Times(ind>0));
                    UVW_ = UVW(ind>0,:)./wavelength;
                   
                    nW_{ch}{conf} = double(Weights_ch(ind>0,j));
                    y_{ch}{conf} = double(I(ind>0) .* nW_{ch}{conf});
                    
		    meanval = 10*mean(abs(y_{ch}{conf}));
                    extra_ind =  (abs(y_{ch}{conf})>0).*(abs(y_{ch}{conf}) <meanval);

                    y_{ch}{conf} =  y_{ch}{conf}(extra_ind>0);
                    numel(y_{ch}{conf})

                    time_{ch}{conf} = (time_{ch}{conf} (extra_ind>0));
                    nW_{ch}{conf} = nW_{ch}{conf} (extra_ind>0);
                    UVW_ = UVW_(extra_ind>0,:);
                else
                    ind = (temp==0).*(abs(I)>0);
                    time_{ch}{conf} = double(Times(ind>0));
                    nW_{ch}{conf} = double(Weights_ch(ind>0,j));
                    y_{ch}{conf} = double(I(ind>0)) .* nW_{ch}{conf};
                    UVW_ = UVW(ind>0,:)./wavelength;
                end
                u_{ch}{conf}  = UVW_(:,1);
                v_{ch}{conf}  = - UVW_(:,2);
                w_{ch}{conf} = UVW_(:,3);UVW_=[];
                maxBaselinePerCh(ch,conf) = max(sqrt(u_{ch}{conf}.^2 + v_{ch}{conf}.^2));
                pos_{ch}(conf) = length(y_{ch}{conf});
                ch = ch + 1;
            end
        end
    end
end
clear u v w y ;
%% Setting imaging params
maxBaseline = max(max(maxBaselinePerCh));
theta = 1 / (2*maxBaseline);
theta_deg = theta * 180 / pi;
theat_arcsec = theta_deg * 3600;

% pixel_size = 0.06  %we write it explicitly to be exact and copy the same value for wsclean
if ~isempty(pixelSize)    &&  pixelSize>0
    SuperResolutionScale = theat_arcsec / pixelSize;
    fprintf('\nINFO: Image resolution: %dx instrumental resolution at the highest channel. \n',SuperResolutionScale);
else
    SuperResolutionScale =2;
    pixelSize = theat_arcsec/SuperResolutionScale ;
    fprintf('\nINFO: Default image resolution: %d asec. \n',pixelSize);
end
ScaleNUFFT =  pi/(maxBaseline * SuperResolutionScale); % convert uv pts in ]-pi,pi[ for NUFFT
nChannels =length(u_);
for i = 1 : nChannels
    for j = 1 : length(u_{i}) % aggregate configurations once done (blocking per config): simple reshape
        u_{i}{j} = u_{i}{j} *ScaleNUFFT;
        v_{i}{j} = v_{i}{j} *ScaleNUFFT;
    end
end
fprintf('\nSaving interleaved data .. \n')
interleaved_channels = split_range_interleaved(Qc, nChannels);
for iSubCube = 1:numel(interleaved_channels)
    iSubCube
    subcube_channels = interleaved_channels{iSubCube};
    u=cell(numel(subcube_channels),1);
    v =cell(numel(subcube_channels),1);
    w=cell(numel(subcube_channels),1);
    time=cell(numel(subcube_channels),1);
    nW=cell(numel(subcube_channels),1);
    y=cell(numel(subcube_channels),1);
    pos =cell(numel(subcube_channels),1);
    for j =1:numel(subcube_channels)
        u{j,1} =u_{subcube_channels(j)};          u_{subcube_channels(j)} =[];
        v{j,1} =v_{subcube_channels(j)};          v_{subcube_channels(j)} =[];
        w{j,1} =w_{subcube_channels(j)};          w_{subcube_channels(j)} =[];
        time{j,1} = time_{subcube_channels(j)};   time_{subcube_channels(j)} =[];
        nW{j,1} = nW_{subcube_channels(j)};       nW_{subcube_channels(j)} =[];
        y{j,1} =y_{subcube_channels(j)};          y_{subcube_channels(j)}=[];
        pos{j,1}=pos_{subcube_channels(j)};     pos_{subcube_channels(j)}=[];
        
    end
    save([file2save,num2str(iSubCube),'.mat'],'-v7.3', 'y', 'u', 'v','w', 'nW', 'time', 'pos', 'pixelSize');
    clear u v w y time nW pos ;
end






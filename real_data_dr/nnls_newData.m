function [norm_res] = nnls_newData(datadir, ch, reduction_version, fouRed_gamma, fouRed_type, wterm, levelG, levelC, lowRes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

addpath ../lib/utils/
addpath ../fouRed/
addpath ../lib/operators
addpath ../lib/nufft

fprintf('Channel number: %d\n', ch)
fprintf('Reduction version %d\n', reduction_version);
if fouRed_type == 1
    typeStr = 'perc';
elseif fouRed_type == 2
    typeStr = 'th';
end
fprintf('Low resolution: %d\n', lowRes);

if lowRes
    Nx = 2560; % 2560;
    Ny = 2560; % 1536;
else
    Nx = 4096; % 2560;
    Ny = 4096; % 1536;
end

%% Config parameters
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 8; % number of neighbours for nufft
Ky = 8; % number of neighbours for nufft

% parameter NNLS
param_nnls.verbose = 2;       % print log or not
param_nnls.rel_obj = 1e-6;    % stopping criterion
param_nnls.max_iter = 1000;     % max number of iterations 1000
param_nnls.sol_steps = [inf]; % saves images at the given iterations
param_nnls.beta = 1;
    
if ~wterm
    [Ap, Atp, ~, ~] = op_nufft([0, 0], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);
end

if wterm
    if lowRes
        DRfilename = [datadir, '/ESO137_LOW_DR_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '=', num2str(ch), '.mat'];
    else
        DRfilename = [datadir, '/ESO137_DR_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '=', num2str(ch), '.mat'];
    end
    fprintf('Read dimensionality reduction file: %s\n', DRfilename)
    tmp = load(DRfilename, 'H', 'W', 'yT', 'T', 'aW', 'Wm', 'A', 'At', 'epsilon');
else
    if lowRes
        DRfilename = [datadir, '/ESO137_LOW_DR_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '=', num2str(ch), '.mat'];
    else
        DRfilename = [datadir, '/ESO137_DR_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '=', num2str(ch), '.mat'];
    end
    fprintf('Read dimensionality reduction file: %s\n', DRfilename)
    tmp = load(DRfilename, 'H', 'W', 'yT', 'T', 'aW', 'Wm', 'epsilon');
end
Hp = tmp.H{1,1};
Wp = tmp.W{1,1};
yTp = tmp.yT{1,1};
Tp = tmp.T{1,1};
aWp = tmp.aW{1,1};
Wmp = tmp.Wm{1,1};
eps_normy = tmp.epsilon{1,1};

if wterm
    Ap = tmp.A;
    Atp = tmp.At;
end
% if reduction_version == 2
%     for j = 1:length(Hp)
%         if usingPrecondition
%             aWp{j} = Tp{j};
%         end
%     end
% end
clear tmp
fprintf('\nDR file of channel number %d has been read\n', ch)
        
[~, norm_res] = fb_dr_nnls(yTp{1}, Ap, Atp, Hp{1}, Wp{1}, Tp{1}, Wmp{1}, param_nnls, reduction_version);
epsilon{1}{1} = norm_res;
fprintf('Effective channel %d, estimated epsilon from NNLS: %f\n', ch, norm_res)
fprintf('Effective channel %d, estimated epsilon from norm(y): %f\n', ch, eps_normy{1})

if wterm
    if lowRes
        DRfilename = [datadir,'/ESO137_LOW_NNLS_DR_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '=', num2str(ch), '.mat'];
    else
        DRfilename = [datadir,'/ESO137_NNLS_DR_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '=', num2str(ch), '.mat'];
    end
    if ~isfile(DRfilename)
        save(DRfilename, '-v7.3', 'epsilon');
    end
else
    if lowRes
        DRfilename = [datadir,'/ESO137_LOW_NNLS_DR_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '=', num2str(ch), '.mat'];
    else
        DRfilename = [datadir,'/ESO137_NNLS_DR_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '=', num2str(ch), '.mat'];
    end
    if ~isfile(DRfilename)
        save(DRfilename, '-v7.3', 'epsilon');
    end
end
end


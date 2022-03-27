function func_compute_residuals(datadir, rsltdir, ch, subInd, reduction_version, realdatablocks, fouRed_gamma, fouRed_type)

if fouRed_type == 1
    typeStr = 'perc';
elseif fouRed_type == 2
    typeStr = 'th';
end

addpath ../fouRed;
addpath ../lib/;
addpath ../lib/operators/;
addpath ../lib/RI-measurement-operator/nufft/;
addpath ../lib/utils/;
addpath ../lib/SARA-dictionary/src;
addpath ../src_mnras/;
addpath ../src_mnras/spmd/;
addpath ../src_mnras/spmd/dr/;
addpath ../src_mnras/spmd/weighted/;

fprintf('Channel number %d\n', ch);
fprintf('Reduction version %d\n', reduction_version);
fprintf('Data blocks: %d\n', realdatablocks);

param_real_data.image_size_Nx = 2560; % 2560;
param_real_data.image_size_Ny = 1536; % 1536;
nChannels = length(ch); % total number of "virtual" channels (i.e., after
% concatenation) for the real dataset considered
% nBlocks = realdatablocks;        % number of data blocks (needs to be known beforehand,
% quite restrictive here), change l.70 accordingly
% klargestpercent = 20;

%% Config parameters
Nx = param_real_data.image_size_Nx;
Ny = param_real_data.image_size_Ny;
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 7; % number of neighbours for nufft
Ky = 7; % number of neighbours for nufft

[A, At, ~, ~] = op_nufft([0, 0], [Ny Nx], [Ky Kx], [oy * Ny ox * Nx], [Ny / 2 Nx / 2]);

% model image
cubefilename = [rsltdir, 'xsol_hs_ddr.fits'];
fprintf('Read model cube file: %s\n', cubefilename);
cube = fitsread(cubefilename);
ci = size(cube, 3);

% nuclear-norm and l21-norm
xhatm = reshape(cube, numel(cube) / ci, ci);
[~, S0, ~] = svd(xhatm, 'econ');
nuclear_norm = norm(diag(S0), 1);
l21_norm = sum(sqrt(sum(abs(xhatm).^2, 2)), 1);
fprintf('\nNuclear-norm = %e, l21-norm = %e\n', nuclear_norm, l21_norm);

% read data and compute residuals
norm_residual_check_c = 0;
norm_residual_check_a = 0;
for i = 1:ci
    % read data
    fprintf('\nChannel number: %d\n', ch(i));
    DRfilename = [datadir, '/CYG_DR_cal_', num2str(realdatablocks), 'b_ind', num2str(subInd(1)), '_', num2str(subInd(end)), '_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '=', num2str(ch(i)), '.mat'];
    fprintf('Read dimensionality reduction file: %s\n', DRfilename);
    tmp = load(DRfilename, 'H', 'W', 'yT', 'T', 'aW', 'Wm');
    H{1, 1} = tmp.H{1, 1};
    W{1, 1} = tmp.W{1, 1};
    yT{1, 1} = tmp.yT{1, 1};
    T{1, 1} = tmp.T{1, 1};
    aW{1, 1} = tmp.aW{1, 1};
    Wm{1, 1} = tmp.Wm{1, 1};
    fprintf('\nDR file of channel number %d has been read\n', ch(i));
    % oversampled size
    R = length(H{1});

    Fx = A(cube(:, :, i));
    for j = 1:R
        r2 = T{1}{j} .* (H{1}{j} * Fx(W{1}{j}));
        % norm of residual
        norm2_res{i}{j} = sum(power(abs(r2 - yT{1}{j}), 2));     % faster than built-in norm
        if (realdatablocks == 2 && j == 1) || (realdatablocks == 9 && (j == 1 || j == 2))
            norm_residual_check_c = norm_residual_check_c + norm2_res{i}{j};
        else
            norm_residual_check_a = norm_residual_check_a + norm2_res{i}{j};
        end
    end
end

norm_residual = sqrt(norm_residual_check_c + norm_residual_check_a);
norm_residual_check_c = sqrt(norm_residual_check_c);
norm_residual_check_a = sqrt(norm_residual_check_a);

fprintf('\nresidual = %e\n', norm_residual);
fprintf('residual_c = %e, residual_a = %e\n', norm_residual_check_c, norm_residual_check_a);

for i = 1:ci
    for j = 1:length(norm2_res{i})
        if realdatablocks == 2 && j == 1
            norm_residual_ic = norm2_res{i}{j};
        else
            norm_residual_ia = norm2_res{i}{j};
        end
    end
    norm_residual_i = sqrt(norm_residual_ic + norm_residual_ia);
    norm_residual_ic = sqrt(norm_residual_ic);
    norm_residual_ia = sqrt(norm_residual_ia);
    fprintf('Channel: %d, residual = %e, residual_c = %e, residual_a = %e\n', ch(i), norm_residual_i, norm_residual_ic, norm_residual_ia);
end
save('residual2.mat', '-v7.3', 'norm2_res');

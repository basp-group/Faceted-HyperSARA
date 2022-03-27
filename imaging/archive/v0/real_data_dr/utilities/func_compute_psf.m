function func_compute_psf(ch, subInd, reduction_version, realdatablocks, fouRed_gamma, fouRed_type)

addpath ../fouRed;
addpath ../lib/;
addpath ../lib/operators/;
addpath ../lib/RI-measurement-operator/nufft/;
addpath ../lib/utils/;
addpath ../lib/faceted-wavelet-transform/src;
addpath ../src_mnras/;
addpath ../src_mnras/spmd/;
addpath ../src_mnras/spmd/dr/;
addpath ../src_mnras/spmd/weighted/;

if fouRed_type == 1
    typeStr = 'perc';
elseif fouRed_type == 2
    typeStr = 'th';
end

fprintf('Channel number %d\n', ch);
fprintf('Index number %d\n', subInd);
fprintf('Reduction version %d\n', reduction_version);
fprintf('Data blocks: %d\n', realdatablocks);
if fouRed_type == 1
    fprintf('Reduction level: remove %f percentile\n', fouRed_gamma);
elseif fouRed_type == 2
    fprintf('Reduction level: keep %f sigma\n', fouRed_gamma);
end

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
mkdir('./psf/');
for i = 1:nChannels
    fprintf('\nChannel number: %d\n', ch(i));
    DRfilename = ['/lustre/home/shared/sc004/dr_2b_result_real_data/CYG_DR_cal_', num2str(realdatablocks), 'b_ind', num2str(subInd(1)), '_', num2str(subInd(end)), '_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '=', num2str(ch(i)), '.mat'];
    fprintf('Read dimensionality reduction file: %s\n', DRfilename);
    tmp = load(DRfilename, 'H', 'W', 'T');
    H = tmp.H{1, 1};
    W = tmp.W{1, 1};
    T = tmp.T{1, 1};
    fprintf('\nDR file of channel number %d has been read\n', ch(i));

    No = size(W{1}, 1);
    dirac2D = zeros(Ny, Nx);
    dirac2D(Ny / 2, Nx / 2) = 1;
    Fx = A(dirac2D);
    g2 = zeros(No, 1);
    for j = 1:length(H)
        r2 = T{j} .* (H{j} * Fx(W{j}));
        g2(W{j}) = g2(W{j}) + H{j}' * (T{j} .* r2);
    end
    psf = real(At(g2));
    psffilename = ['./psf/psf_', num2str(realdatablocks), 'b_ind', num2str(subInd(1)), '_', num2str(subInd(end)), '=', num2str(ch(i)), '.fits'];
    if ~isfile(psffilename)
        fitswrite(psf, psffilename);
    end
    fprintf('\nPSF of channel number %d has been computed\n', ch(i));
end
end

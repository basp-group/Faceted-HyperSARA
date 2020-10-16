function uv_bin_sampling(inputdir, outputdir, nb_aggch, chInd, subProb, wterm, levelG, levelC, lowRes)

addpath ../lib/utils/
addpath ../fouRed/
addpath ../lib/operators
addpath ../lib/nufft

addpath ../../AD_build_wproj_op

speed = 299792458;  % light speed

fprintf('Channel number: %d\n', chInd);
fprintf('Subprob number: %d\n', subProb);
fprintf('W term: %d\n', wterm);
if wterm
    fprintf('Energy and chirp level: %f, %f\n', levelG, levelC);
end
fprintf('Low resolution: %d\n', lowRes);

% load data
datafilename = [inputdir,'ESO137_FULL.mat'];
fprintf("Read file: %s\n", datafilename)
load(datafilename,'weights_ch','vis','Freqs','uvw');
% load final flag
flagfilename = [inputdir,'final_flag.mat'];
fprintf("Read file: %s\n", flagfilename)
load(flagfilename, 'bmax', 'dl', 'final_flag')

%% Image size
if lowRes
    Nx = 2560;
    Ny = 2560;
else
    Nx = 4096;
    Ny = 4096;
end

if lowRes
    dl = dl/1.6;
end
fprintf("dl: %f\n", dl)

if lowRes
    imPixelSize = 4096*1/max(Nx, Ny);    % resolution in arcsecond
else
    imPixelSize = 1;    % resolution in arcsecond
end
fprintf("Cell size: %f arc second\n", imPixelSize)

% Config parameters
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
if wterm
    Kx = 7; % number of neighbours for nufft
    Ky = 7; % number of neighbours for nufft
else
    Kx = 7; % number of neighbours for nufft
    Ky = 7; % number of neighbours for nufft
end

% parameter NNLS
param_nnls.verbose = 2;       % print log or not
param_nnls.rel_obj = 1e-5;    % stopping criterion
param_nnls.max_iter = 500;     % max number of iterations 1000
param_nnls.sol_steps = [inf]; % saves images at the given iterations
param_nnls.beta = 1;

% Effective channel configuration (data-related)
bands = [33:960];
if wterm
    numworkers = 36;
    delete(gcp('nocreate'))
    cirrus_cluster = parcluster('local');   % slurm
    parpool(cirrus_cluster, numworkers, 'IdleTimeout', Inf);
end
for l = 1:length(chInd) %1:nb_effch
    
    j = chInd(l);
    
%     yw = [];
%     nW = [];
%     uvw_rad = [];
    nb_points = 0;
    for k = 1:length(subProb)
        subk = subProb(k);
        if j==0
            ind = subk;
        else        
            ind = (j-1)*nb_aggch+subk+bands(1)-1;
        end
        fprintf("Effective channel: %d, channel:%d\n", j, ind)
        % uvw
        uvw1 = uvw(~final_flag{ind},:);
        uvw1(:,2) = -uvw1(:,2);   % change to coordinates [u,-v,w]
        wave = speed/Freqs(ind);
%         uvw_rad = [uvw_rad; uvw1./wave*pi/(bmax*dl)];
%         nW = [nW; double(weights_ch(~final_flag{ind})')];
%         yw = [yw; double(vis(~final_flag{ind},ind))];
        uvw_rad = uvw1./wave*pi/(bmax*dl);
        nW = double(weights_ch(~final_flag{ind})');
        yw = double(vis(~final_flag{ind},ind));
        if wterm
            [G,opA,opAt,Lnorm,levelG,supports]= get_wplane_wprojection_G(uvw_rad*bmax*dl/pi, Nx, Ny, imPixelSize, levelG, levelC, nW, bmax, dl);
            A = opA{1};
            At = opAt{1};
        else
            [A, At, G, ~] = op_nufft([uvw_rad(:,2) uvw_rad(:,1)], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);   % nufft using uvw in rad
            G = spdiags(nW, 0, size(uvw_rad,1), size(uvw_rad,1)) * G;
        end
        Wl = any(G, 1).';
        G = G(:, Wl);
        Gw{l}{k} = G;
        W{l}{k} = Wl;
        yT{l}{k} = nW .* yw;
        nb_points = length(yw) + nb_points;
        clear G Wl yw nW
    end
%     fprintf("nb of points:%d\n", length(yw(:)))
    fprintf("nb of points:%d\n", nb_points)
    % vis
%     yw = nW .* yw;    % natural-weighted vis
    
%     if length(yw(:)) ~= length(nW(:)) || length(yw(:)) ~= size(uvw_rad,1) || length(nW(:)) ~= size(uvw_rad,1)
%         fprintf("Error: Dimension not consistent!")
%     end
%     if wterm
%         [G,opA,opAt,Lnorm,levelG,supports]= get_wplane_wprojection_G(uvw_rad*bmax*dl/pi, Nx, Ny, imPixelSize, levelG, levelC, nW, bmax, dl);
%         A = opA{1};
%         At = opAt{1};
%     else
%         [A, At, G, ~] = op_nufft([uvw_rad(:,2) uvw_rad(:,1)], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);   % nufft using uvw in rad
%         G = spdiags(nW, 0, size(uvw_rad,1), size(uvw_rad,1)) * G;
%     end
    
%     Wl = any(G, 1).';
%     G = G(:, Wl);
    
%     if wterm
%         fitswrite(uvw_rad, ['uv_bin', num2str(nb_aggch), '_sampling_ch', num2str(j), '.fits'])
%     end
%     clear uvw_rad nW
    
%     eps_normy = 0.01 * norm(yw(:));
%     fprintf("Epsilon estimated from norm(y): %e\n", eps_normy)
%     [~, norm_res] = fb_nnls_blocks(yw, A, At, G, Wl, param_nnls);
%     fprintf("Epsilon estimated from NNLS: %e\n", norm_res)

    
    epsilon{l} = cell(length(subProb),1);
    epsilon_normy{l} = cell(length(subProb),1);
    for k = 1:length(subProb)
        fprintf('solving for band %i\n\n', k)
        epsilon_normy{l}{k} = 0.01 * norm(yT{l}{k});
        fprintf("Epsilon estimated from norm(y): %e\n", epsilon_normy{l}{k})
        [~,epsilon{l}{k}] = fb_nnls_blocks(yT{l}{k}, A, At, Gw{l}{k}, W{l}{k}, param_nnls);
        fprintf("Epsilon estimated from NNLS: %e\n", epsilon{l}{k})
    end

    
%     Gw{1}{1} = G;
%     W{1}{1} = Wl;
%     yT{1}{1} = yw;
%     epsilon{1}{1} = norm_res;
%     epsilon_normy{1}{1} = eps_normy;
    
    if wterm
        if lowRes
            filename = [outputdir,'/ESO137_LOW_w', num2str(levelG), '_', num2str(levelC), '=', num2str(j), '.mat'];
        else
            filename = [outputdir,'/ESO137_w', num2str(levelG), '_', num2str(levelC), '=', num2str(j), '.mat'];
        end
        if ~isfile(filename)
            save(filename, '-v7.3', 'A','At','Gw', 'W', 'yT', 'epsilon', 'epsilon_normy');
        end
    else
        if lowRes
            filename = [outputdir,'/ESO137_LOW=', num2str(j), '.mat'];
        else
            filename = [outputdir,'/ESO137=', num2str(j), '.mat'];
        end
        if ~isfile(filename)
            save(filename, '-v7.3', 'Gw', 'W', 'yT', 'epsilon', 'epsilon_normy');
        end
    end
    
end

fprintf("\nUV samplings are finished binning!\n")
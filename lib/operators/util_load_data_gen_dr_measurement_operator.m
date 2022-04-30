function [A, At, H, W, aW, Lambda, gridded_y, gridded_noise] = util_load_data_gen_dr_measurement_operator(y, dataFilename, effChans2Image, ...
                                                                                                          nDataSets, Nx, Ny, param_nufft, param_wproj, preproc_noise_sigma, preproc_filenames)
    % Build the measurement operator incorporating visibility gridding
    % for a given uv-coverage at pre-defined frequencies.
    %
    % Parameters
    % ----------
    % y : cell of cell of complex[:]
    %    Cell containing the data vectors (Stokes I) to be gridded.
    % dataFilename: function handle
    %     Filenames of the data files to load ``u``, ``v`` and ``w`` coordinates and
    %     ``nW`` the weights involved in natural weighting.
    % effChans2Image: cell
    %     Indices of the  physical (input) channels, combined to image the
    %     effective (ouput) channel.
    % nDataSets :int
    %     Number of datasets per physical (input) channel.
    % Nx : int
    %     Image dimension (x-axis).
    % Ny : int
    %     Image dimension (y-axis).
    % param_nufft : struct
    %     Structure to configure NUFFT.
    % param_wproj : struct
    %     Structure to configure `w`-projection.
    % param_precond : struct
    %     Structure to configure the preconditioning matrices.
    % preproc_noise_sigma: cell of doubles
    %     Standard deviation values of the noise in the physical (input)
    %     channels.
    % preproc_filenames: struct
    %     Struct composed of two fields: (i) ``preproc_filenames.dde``: anonymous
    %     functions to get DDE calibration kernel files, expected
    %     variable to be loaded ``DDEs``. (ii) ``preproc_filenames.noise_std``:
    %     anonymous functions to get the noise-related files,  expected
    %     variable to be loaded (optional) ``RESIDUAL``  residual
    %     visibility vector from a preprocessing step, assumed to be a noise
    %     vector.
    %
    % Returns
    % -------
    % A : function handle
    %     Function to compute the rescaled 2D Fourier transform involved
    %     in the emasurement operator.
    % At : function handle
    %     Function to compute the adjoint of ``A``.
    % H : cell of cell of complex[:]
    %     Cell containing the holographic matrix for each channel, and
    %     each data block within a channel. H is a lower-triangular
    %     matrix to optimise memory requirements.
    % W : cell of cell of double[:]
    %     Cell containing the selection vector for each channel, and
    %     data block within a channel.
    % aWw : cell of cell of double[:]
    %     Cell containing the preconditioning vectors for each channel, and
    %     data block within a channel.
    % Lambda: cell of cell of double
    %     Cell containing the diagonal weighting matrix involved in data
    %     dimensionality reduction via visibility gridding.
    % gridded_y: cell of cell of complex[:]
    %     Cell containing the reduced data vectors which corresponds to the gridded visibilities.
    % gridded_noise: struct
    %     struct containing estimated l2bounds and  noise std, from
    %     gridding a noise data vector.

    %%
    speed_of_light = 299792458;
    % check if data pre-processing exist
    % ddes
    if ~exist('preproc_filenames', 'var')
        ddesfilename = [];
    elseif ~isfield(preproc_filenames, 'dde')
        ddesfilename = [];
    else
        ddesfilename = preproc_filenames.dde;
    end
    % residual data
    if ~exist('preproc_filenames', 'var')
        filename_noise_std = [];
    elseif ~isfield(preproc_filenames, 'noise_std')
        filename_noise_std = [];
    else
        filename_noise_std = preproc_filenames.noise_std;
    end
    % noise estimates
    if  ~exist('preproc_dr_noise_sigma', 'var')
        preproc_noise_sigma = [];
    end

    % dims.
    param_nufft.N = [Ny Nx];
    param_nufft.Nn = [param_nufft.Ky param_nufft.Kx];
    param_nufft.No = [param_nufft.oy * Ny param_nufft.ox * Nx];
    param_nufft.Ns = [Ny / 2  Nx / 2];

    % init.
    nchans = numel(effChans2Image);
    H = cell(nchans, 1);
    W = cell(nchans, 1);
    aW = cell(nchans, 1);
    Lambda = cell(nchans, 1);
    gridded_y = cell(nchans, 1);
    gridded_noise = [];

    % get Fourier operators
    [A, At] = op_p_nufft_wproj_dde(param_nufft);
    flag_apply_imaging_weights = param_nufft.flag_apply_imaging_weights;
    % get holographic matrix and other vars
    for ifc = 1:nchans
        % INFO!!  no -precond for DR  and 1 single block per channel for now
        aW{ifc}{1} = 1; % no precond ...
        gridded_noise_vect =  sparse(prod(param_nufft.No), 1);
        % measurement operator initialization
        H{ifc}{1} = sparse(prod(param_nufft.No), prod(param_nufft.No));
        gridded_y{ifc}{1} = sparse(prod(param_nufft.No), 1);
        for iCh = 1:numel(effChans2Image{ifc})
            fprintf('Channel %d. ', effChans2Image{ifc}(iCh));
            for idSet = 1:nDataSets
                % fprintf('Channel %d, dataset id %d.  ',effChans2Image{ifc}(iCh),idSet)
                if flag_apply_imaging_weights
                    load(dataFilename(idSet, effChans2Image{ifc}(iCh)), 'u', 'v', 'w', 'nW', 'frequency', 'nWimag');
                    try
                        nW = double(double(nW) .* sqrt(double(nWimag)));
                        nWimag = [];
                    catch
                        flag_apply_imaging_weights = false;
                    end
                else
                    load(dataFilename(idSet, effChans2Image{ifc}(iCh)), 'u', 'v', 'w', 'nW', 'frequency');
                end
                % u v  are in units of the wavelength and will be normalised between [-pi,pi] for the NUFFT
                u = double(u(:)) * pi / param_wproj.halfSpatialBandwidth;
                v = -double(v(:)) * pi / param_wproj.halfSpatialBandwidth;
                w = double(w(:));
                nW = double(nW(:));
                DDEs = [];
                if ~isempty(ddesfilename)
                    tmpfile = ddesfilename(idSet, effChans2Image{ifc}(iCh));
                    if isfile(tmpfile) % calibration kernels
                        cal_solutions = load(tmpfile);
                        if isfield(cal_solutions, 'DIEs')
                            if ~isequal(size(nW), size(cal_solutions.DIEs))
                                try
                                    nW = nW .* cal_solutions.DIEs(abs(cal_solutions.DIEs) ~= 1);
                                catch
                                    fprintf('SEVERE: DIEs dimensions are not compatible, no DIEs applied.\n');
                                end
                            else
                                nW  = nW .* cal_solutions.DIEs; % include DIEs directly with the weights
                            end
                        end
                        if isfield(cal_solutions, 'DDEs')
                            DDEs = cal_solutions.DDEs;
                        end
                        clear cal_solutions;
                    else
                        fprintf('\nWARNING: Calibration solutions not found: %s ', tmpfile);
                    end
                end
                %% get the de-gridding matrix
                if isempty(DDEs)
                    [~, ~, G, W_curr] = op_p_nufft_wproj_dde(param_nufft, [{v} {u}], {w}, {nW}, param_wproj);
                else
                    [~, ~, G, W_curr] = op_p_nufft_wproj_dde(param_nufft, [{v} {u}], {w}, {nW}, param_wproj, {DDEs});
                end
                clear u v w DDEs;

                G_curr = sparse(prod(param_nufft.No), size(G{1}, 1));
                G_curr(W_curr{1}, :) = G{1}';
                clear G;
                %% grid data
                gridded_y{ifc}{1} = sparse(gridded_y{ifc}{1}  + G_curr * (y{ifc}{iCh}{idSet} .* (abs(nW) > 0)));

                %% grid noise vector
                try % existing noise vector
                    % load & grid noise vector
                    tmpfile = filename_noise_std(idSet, effChans2Image{ifc}(iCh));
                    noise_vect = load(tmpfile, 'RESIDUAL');
                    gridded_noise_vect = sparse(gridded_noise_vect + G_curr * noise_vect.RESIDUAL);
                    clear noise_vect;
                catch % generate random noise vect
                    if ~isempty(preproc_noise_sigma) % estimated noise std
                        if ~isempty(preproc_noise_sigma{ifc}{iCh}{idSet})
                            sigma_cont_vis = preproc_noise_sigma{ifc}{iCh}{idSet};
                        end
                    else % theoretical noise std
                        sigma_cont_vis = 1;
                    end
                    % generate & grid noise vector
                    nmeas_blk = numel(y{ifc}{iCh}{idSet});
                    gridded_noise_vect = sparse(gridded_noise_vect + ...
                                                G_curr * (sigma_cont_vis .* (randn(nmeas_blk, 1) + 1i * randn(nmeas_blk, 1)) / sqrt(2)));
                end

                %% update H
                H{ifc}{1} = H{ifc}{1} +  tril(G_curr * G_curr'); % keep lower half of H only for memory reasons
                clear G_curr;
            end
        end
        %% get diag elements of H
        diagH = diag(H{ifc}{1});

        %% get active Fourier modes
        W{ifc}{1} = (diagH > 0);

        %% update Lambda
        Lambda{ifc}{1} = 1 ./ sqrt(diagH(W{ifc}{1}));
        clear diagH;
        %% apply Lambda to data
        gridded_y{ifc}{1} = gridded_y{ifc}{1}(W{ifc}{1}) .* Lambda{ifc}{1};

        %% get noise stats when available residual data
        gridded_noise_vect =  gridded_noise_vect(W{ifc}{1}) .*  Lambda{ifc}{1}; % apply Lambda to residual data
        gridded_noise.l2bounds{ifc}{1} = sqrt(sum(abs(gridded_noise_vect).^2));
        gridded_noise.sigma{ifc} = std(full(gridded_noise_vect));
        clear gridded_noise_vect;

        %% Update H and its diag elements.
        H{ifc}{1} = H{ifc}{1}(W{ifc}{1}, W{ifc}{1});
        if istril(H{ifc}{1})
            % Update its diag elements to apply it & its transpose conj. in the M.O.
            for idiag = 1:size(H{ifc}{1}, 1)
                H{ifc}{1}(idiag, idiag) = 0.5 * H{ifc}{1}(idiag, idiag);
            end
        end
        fprintf('\nH is built succesfully.');
        whos H;
    end

end

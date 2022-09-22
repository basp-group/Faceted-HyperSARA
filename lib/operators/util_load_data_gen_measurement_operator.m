function [A, At, G, W, aW] = util_load_data_gen_measurement_operator(dataFilename, effChans2Image, ...
                                                                     nDataSets, Nx, Ny, param_nufft, param_wproj, param_precond, ddesfilename)
    % Build the measurement operator for a given uv-coverage at pre-defined
    % frequencies.
    %
    % Parameters
    % ----------
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
    %     Structure to configure w-projection.
    % param_precond : struct
    %     Structure to configure the preconditioning matrices.
    % ddesfilename: function handle
    %     Filenames of the DDE calibration kernels in the spatial Fourier domain
    %     from a pre-processing step to be incorporated in the measurement
    %     operator. Expected variable to be loaded ``DDEs``.
    %

    % Returns
    % -------
    % A : function handle
    %     Function to compute the rescaled 2D Fourier transform involved
    %     in the emasurement operator.
    % At : function handle
    %     Function to compute the adjoint of ``A``.
    % G : cell of cell of complex[:]
    %     Cell containing the trimmed-down interpolation kernels for each
    %     channel, and each data block within a channel.
    % W : cell of cell of double[:]
    %     Cell containing the selection vector for each channel, and
    %     data block within a channel.
    % aW : cell of cell of double[:]
    %     Cell containing the preconditioning vectors for each channel, and
    %     data block within a channel.
    %%
    speed_of_light = 299792458;

    if ~exist('ddesfilename', 'var')
        ddesfilename = [];
    end
    param_nufft.N = [Ny Nx];
    param_nufft.Nn = [param_nufft.Ky param_nufft.Kx];
    param_nufft.No = [param_nufft.oy * Ny param_nufft.ox * Nx];
    param_nufft.Ns = [Ny / 2 Nx / 2];

    nchans = numel(effChans2Image);
    G = cell(nchans, 1);
    W = cell(nchans, 1);
    aW = cell(nchans, 1);

    % get fft operators
    [A, At, ~, ~] = op_p_nufft_wproj_dde(param_nufft);
    flag_apply_imaging_weights = param_nufft.flag_apply_imaging_weights;

    for ifc = 1:nchans
        nBlocks  = numel(effChans2Image{ifc}) * nDataSets;
        % set the blocks structure
        aW{ifc} = cell(nBlocks, 1);
        G{ifc} = cell(nBlocks, 1);
        W{ifc} = cell(nBlocks, 1);
        ct = 1;
        for iCh = 1:numel(effChans2Image{ifc})
            for idSet = 1:nDataSets
                if flag_apply_imaging_weights
                    load(dataFilename(idSet, effChans2Image{ifc}(iCh)), 'u', 'v', 'w', 'nW', 'frequency', 'nWimag');
                    try
                        nW = double(double(nW) .* (double(nWimag)));
                        nWimag = [];
                    catch
                        flag_apply_imaging_weights = false;
                    end
                else
                    load(dataFilename(idSet, effChans2Image{ifc}(iCh)), 'u', 'v', 'w', 'nW', 'frequency');
                end
                wavelength = speed_of_light / frequency;
                % u v  are in units of the wavelength and will be normalised between [-pi,pi] for the NUFFT
                u = double(u(:)) * pi / param_wproj.halfSpatialBandwidth;
                v = -double(v(:)) * pi / param_wproj.halfSpatialBandwidth;
                w = double(w(:));
                nW = double(nW(:));
                DDEs = [];
                if ~isempty(ddesfilename)
                    tmpfile = ddesfilename(idSet, effChans2Image{ifc}(iCh));
                    if isfile(tmpfile) % load calibration kernels
                        cal_solutions = load(tmpfile);
                        if isfield(cal_solutions, 'DIEs')
                            if ~isequal(size(nW), size(cal_solutions.DIEs))
                                fprintf('SEVERE: DIEs Dimensions are not compatible, no DIEs applied.\n');
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

                % compute uniform weights (sampling density) for the preconditioning
                aW{ifc}{ct, 1} = util_gen_preconditioning_matrix(u, v, param_precond);
                % measurement operator initialization
                if isempty(DDEs)
                    [~, ~, G_curr, W_curr] = op_p_nufft_wproj_dde(param_nufft, [{v} {u}], {w}, {nW}, param_wproj);
                else
                    [~, ~, G_curr, W_curr] = op_p_nufft_wproj_dde(param_nufft, [{v} {u}], {w}, {nW}, param_wproj, {DDEs});
                end
                clear v u w nW DDEs;
                G{ifc}{ct, 1} = G_curr{1};
                clear G_curr;
                W{ifc}{ct, 1} = W_curr{1};
                clear W_curr;
                ct = ct + 1;
            end

        end
    end

end

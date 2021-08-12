function [A, At, G, W] = op_p_nufft_wproj_dde(p,w, nufft,wproj, param)

% Create the nonuniform gridding matrix and fft operators to be used for
% parallel processing
%
% in:
% p{:}[2] - nonuniformly distributed frequency location points for each
%           cell member which will be treated in parallel
% N[2]    - size of the reconstruction image
% Nn[2]   - size of the kernels (number of neighbors considered on each direction)
% No[2]   - oversampled fft from which to recover the non uniform fft via
%           kernel convolution
% Ns[2]   - fft shift
%
% out:
% A[@]          - function handle for direct operator
% At[@]         - function handle for adjoint operator
% G{:}[:][:]    - convolution kernel matrix (small) associated with each
%               patch in the fourier plane
% W{:}          - mask of the values that contribute to the convolution
% Gw[:][:]      - global convolution kernel matrix
%%
%A. Onose, A. Dabbech, Y. Wiaux - An accelerated splitting algorithm for radio-interferometric %imaging: when natural and uniform weighting meet, MNRAS 2017, arXiv:1701.01748
%https://github.com/basp-group/SARA-PPD
%%
N =nufft.N;
Nn=nufft.Nn;
No=nufft.No;
Ns=nufft.Ns;




if ~exist('param', 'var')
    param = struct();
end
if ~isfield(param, 'use_nufft_blocks'), param.use_nufft_blocks = 1; end
if ~isfield(param, 'gen_only_fft_op'), param.gen_only_fft_op = 0; end


R = size(p, 1);
if param.gen_only_fft_op
    [A, At, ~, ~] = op_nufft([0, 0], nufft.N, nufft.Nn, nufft.No, nufft.Ns);
    G = [];
    W = [];
    Gw = [];
else
    
    if ~isfield(nufft, 'nW')
        nufft.nW = cell(length(p), 1);
        for q=1:length(p)
            nufft.nW{q} = ones(length(p{q}(:, 1)), 1);
        end
    end
    if ~param.use_nufft_blocks
        %% compute the overall gridding matrix and its associated kernels
        [A, At, Gw, ~] = op_nufft(cell2mat(p), nufft.N, nufft.Nn, nufft.No, nufft.Ns);

        %% compute small gridding matrices associated with each parallel block
        G = cell(R, 1);
        W = cell(R, 1);

        % block start position
        fprintf('\nComputing block matrices ...\n');
        b_st = 1;
        for q = 1:R
            tstart = tic;
            % current block length
            % the matrix Gw is structured identical to the structure of p thus we 
            % grab it block by block
            b_l = length(p{q});

            % get a block out of the large G and trim it
            Gw(b_st:b_st+b_l-1, :) = spdiags(nufft.nW{q}, 0, b_l, b_l) * Gw(b_st:b_st+b_l-1, :);
            Gw = Gw(b_st:b_st+b_l-1, :);

            %% now trim the zero rows and store a mask in W

            % preallocate W for speed
%             W{q} = false(No(1)*No(2), 1);

            % use the absolute values to speed up the search
%             Gb_a = abs(Gb);

            % check if eack line is entirely zero
%             W{q} = Gb_a' * ones(b_l, 1) ~= 0;
            W{q} = any(Gw, 1).';
            
            % store only what we need from G
            G{q} = Gw(:, W{q});
%             G{q} = Gb;

            % iterate among the blocks
            b_st = b_st+b_l;
            tend = toc(tstart);
            fprintf('Block matrix %d: %ds \n', q, ceil(tend));
        end
    else

        %% compute small gridding matrices associated with each parallel block
        
        %Gw = spalloc(length(cell2mat(p)), No(1)*No(2), 16 * length(cell2mat(p)));
        G = cell(R, 1);
        W = cell(R, 1);

        b_st = 1;
        % block start position
        fprintf('\nComputing block matrices ...\n');
        for q = 1:R

            tstart = tic;
            b_l = length(p{q});


            %% compute the small gridding matrix and its associated kernels
            [~, ~, Gw, ~] = op_nufft([p{q, 1} p{q, 2}], nufft.N, nufft.Nn, nufft.No, nufft.Ns);

            %% now trim the zero rows and store a mask in W

            % preallocate W for speed
            W{q} = false(No(1)*No(2), 1);
            Gw = spdiags(nufft.nW{q}, 0, b_l, b_l) * Gw;
            % check if w correction is needed  
            effBandwidthWterm = max(abs(max(wproj.FoVy, wproj.FoVx).*abs(w{q})));
            wproj.needed =0;
            if effBandwidthWterm > 3*max(wproj.uGridSize,wproj.vGridSize) % hard coded limit of the w bandwidth
                fprintf('\nINFO:W-correction might be needed ..\n')
                wproj.needed =1;
            end
            if (isempty(wproj.do) || wproj.do)  && wproj.needed  
                fprintf('\nINFO: W-correction will be performed\n')
                wproj.do = 1;
              
            else, fprintf('\nINFO: NO W-correction\n')
                  wproj.do =0;
            end
            % perform w-correction
            if  wproj.do
                wproj.paddFourierFactor =nufft.No./nufft.N;
                wproj.gImDims  = nufft.N;
                wproj.supportK=prod(nufft.Nn);
                Gw = Gw.';
                Gw = fun_wprojection(Gw,w{q},[],wproj);
                Gw = Gw';
            end
            
            % check if eack line is entirely zero
            W{q} = any(abs(Gw), 1).';

            % store only what we need from G
            G{q} = Gw(:, W{q});

            %% fill the whole Gw
            %Gw(b_st:b_st+b_l-1, :) = Gb;

            b_st = b_st+b_l;
            tend = toc(tstart);
            fprintf('Block matrix %d: %ds \n', q, ceil(tend));
        end

        [A, At, ~, ~] = op_nufft([0, 0], nufft.N, nufft.Nn, nufft.No, nufft.Ns);
    end
end


end


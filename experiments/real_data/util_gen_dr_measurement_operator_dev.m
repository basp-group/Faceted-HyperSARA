function [A, At, H, W, aW,Sigma,data,noise] = util_gen_dr_measurement_operator_dev(y,u, v,w,nWw, ...
    nchans, Nx, Ny, param_nufft,param_wproj,param_precond, param_blocking,preproc_dr_residuals,ddes)
% Build the measurement operator for a given uv-coverage at pre-defined
% frequencies.
%
% Build the measurement operator for a given uv-coverage at pre-defined
% frequencies.
%
% Parameters
% ----------
% u : array (vector)
%     `u` coordinate.
% v : array (vector)
%     `v` coordinate.
% w : array (vector)
%     `w` coordinate.
% param_precond : struct
%     Structure to configure the preconditioning matrices.
% param_blocking : struct
%     Structure to configure data blocking.
% nchans : int
%     Number of channels.
% Nx : int
%     Image dimension (x-axis).
% Ny : int
%     Image dimension (y-axis).
% Kx : int
%     Dimension interpolation kernel (x-axis).
% Ky : int
%     Dimension interpolation kernel (y-axis).
% ox : int
%     Fourier oversampling factor(x-axis).
% oy : int
%     Fourier oversampling factor(y-axis).
% kernel : string
%     Type of interpolation kernel selected ('kaiser' or 'minmax:tuned').
%
% Returns
% -------
% A : lambda function
%     Lambda function to compute the rescaled 2D Fourier transform involved
%     in the emasurement operator.
% At : lambda function
%     Lambda function to compute the adjoint of ``A``.
% G : cell
%     Cell containing the trimmed-down interpolation kernels for each
%     channel, and each data block within a channel.
% W : cell
%     Cell containing the selection vector for each channel, and
%     data block within a channel.
% aWw : cell
%     Cell containing the preconditioning vectors for each channel, and
%     data block within a channel.
%%
% check if data pre-processing exist
%ddes
if nargin<14
    ddes = [];
end
%residual data
if nargin<13
    ddes = [];
    preproc_dr_residuals =[];
end

% dims.
param_nufft.N =[Ny Nx];
param_nufft.Nn=[param_nufft.Ky param_nufft.Kx];
param_nufft.No= [param_nufft.oy*Ny param_nufft.ox*Nx];
param_nufft.Ns = [Ny/2 Nx/2];

% init.
H = cell(nchans, 1);
W = cell(nchans, 1);
aW = cell(nchans, 1);
Sigma = cell(nchans, 1);
data = cell(nchans, 1);
noise = [];

% get Fourier operators
[A, At] = op_p_nufft_wproj_dde([],1,param_nufft,param_wproj);

% get holographic matrix and other vars
for i = 1:nchans
    %INFO!!  no -precond for DR  and 1 single block per channel for now
    aW{i}{1} = 1;% no precond ...
    
    nW =  nWw{i};
    res_data =  [];
    
    % measurement operator initialization
    H{i}{1} = sparse(prod(param_nufft.No),prod(param_nufft.No));
    
    for j =1:numel(u{i})
        fprintf('\nData Block %d \n',j);
        param_nufft.nW = nW(j);
        if isempty(ddes)
            [~,~, G, W_curr] = op_p_nufft_wproj_dde([v{i}(j) u{i}(j)],w{i}(j),param_nufft,param_wproj);
        else
            [~,~, G, W_curr] = op_p_nufft_wproj_dde([v{i}(j) u{i}(j)],w{i}(j),param_nufft,param_wproj,ddes{i}(j) );
        end
        G_curr = sparse(prod(param_nufft.No),size(G{1},1));
        G_curr(W_curr{1},:)= G{1}'; %#ok<SPRIX>
        clear G;
        
        %% grid residual data when available
        if ~isempty(preproc_dr_residuals)
            if j==1
                res_data = G_curr * preproc_dr_residuals{i}{j};
            else
                res_data =  res_data + G_curr * preproc_dr_residuals{i}{j};
            end
        end
        
        %% grid data 
        if j==1
            data{i}{1} = G_curr * y{i}{j};
        else
            data{i}{1} =  data{i}{1}  + G_curr * y{i}{j};
        end
        
        %% update H
        H{i}{1} = H{i}{1} + tril(G_curr * G_curr'); % keep lower half of H only for memory reasons
        clear Gcurr;
    end
   
    %% get diag elements of H
    diagH = diag(H{i}{1});
    
    %% get active Fourier modes
    W{i}{1} = (diagH>0); 
    
    %% update Sigma
    Sigma{i}{1} = 1./sqrt(diagH(W{i}{1}));
    
    %% apply Sigma to data
    data{i}{1} = data{i}{1}( W{i}{1} ) .* Sigma{i}{1} ; 
    
    %% get noise stats when available residual data 
    if ~isempty(res_data)
        res_data( W{i}{1} ) =  res_data( W{i}{1} ) .*  Sigma{i}{1}; % apply Sigma to residual data
        noise.l2bounds{i}{1} = norm(res_data( W{i}{1} ));
        noise.sigma{i} = std(nonzeros([real(res_data) ; imag(res_data)]));
    end, clear res_data;
   
    
    %% Update H and its diag elements.
    H{i}{1} = H{i}{1}( W{i}{1}, W{i}{1} );
    if istril(H{i}{1})
        % Update its diag elements to apply it & its transpose conj. in the M.O.
        for idiag =1: size(H{i}{1},1)
            H{i}{1}(idiag,idiag) = 0.5* H{i}{1}(idiag,idiag);
        end
    end
    
end

end

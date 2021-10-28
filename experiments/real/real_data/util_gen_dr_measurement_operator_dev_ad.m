function [A, At, H, W, aW,Sigma,data] = util_gen_dr_measurement_operator_dev_ad(y,u, v,w,nWw, ...
    param_precond, param_blocking, nchans, Nx, Ny, param_nufft,param_wproj,param_preproc)
% Build the measurement operator for a given uv-coverage at pre-defined
% frequencies.
%
% Build the measurement operator for a given uv-coverage at pre-defined
% frequencies.
%
% Parameters
% ----------
% u : array (vector)
%     `u` coverage.
% v : array (vector)
%     `v` coverage.
% param_precond : struct
%     Structure to configure the preconditioning matrices.
% param_blocking : struct
%     Structure to configure data blocking.
% fc : array (vector)
%     Frequencies at which the operator needs to be defined.
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
param_nufft.N =[Ny Nx];
param_nufft.Nn=[param_nufft.Ky param_nufft.Kx];
param_nufft.No= [param_nufft.oy*Ny param_nufft.ox*Nx];
param_nufft.Ns = [Ny/2 Nx/2];

H = cell(nchans, 1);
W = cell(nchans, 1);
aW = cell(nchans, 1);
data =  cell(nchans, 1);
Sigma =  cell(nchans, 1);
param.gen_only_fft_op=1;
[A, At, ~,~] = op_p_nufft_wproj_dde([0 0],1,param_nufft,param_wproj,param );

for i = 1:nchans

    % set the blocks structure
    if ~isempty(param_blocking)
        % compute uniform weights (sampling density) for the preconditioning
        aWw = util_gen_preconditioning_matrix(u{i}, v{i}, param_precond);
        % set the weighting matrix, these are the natural weights for real data
        [u{i}, v{i}, ~, ~, aW{i}, nW] = util_gen_block_structure(u{i}, v{i}, aWw{i}, nWw{i}, param_blocking);
        
    else
        for j =1:numel(u{i})
            aW{i}{j,1} = 1;%util_gen_preconditioning_matrix(u{i}{j}, v{i}{j}, param_precond);
        end
        nW =  nWw{i};
    end
    % measurement operator initialization
    % [A, At, G{i}, W{i}] = op_p_nufft_irt([v1 u1], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW, kernel);
   
     %[Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW %(p, N, Nn, No, Ns, ww, param)
    if ~(param_preproc.done && param_preproc.read_G)
        
        H{i}{1} = sparse(prod(param_nufft.No),prod(param_nufft.No));
        W_ = sparse(prod(param_nufft.No),1);
        for j =1:numel(u{i})
            fprintf('\nData Block %d \n',j);
            param_nufft.nW = nW(j);
            [~,~, G, Wcurr] = op_p_nufft_wproj_dde([v{i}(j) u{i}(j)],w{i}(j),param_nufft,param_wproj );
            Gcurr = sparse(prod(param_nufft.No),size(G{1},1));
            Gcurr(Wcurr{1},:)= G{1}'; G =[];
            % gridded data
            tic
            if j==1, data{i}{1} = Gcurr*y{i}{j};
            else, data{i}{1} =  data{i}{1}  + Gcurr*y{i}{j};
            end
            toc
            %H
            tic
            H{i}{1} = H{i}{1} +tril(Gcurr*Gcurr'); % keep lower half only of H
            toc
            
            Gcurr =[];
            W_ = sparse((W_ +Wcurr{1})>0);
            
        end,clear Gcurr;
        
        H{i}{1} = H{i}{1}(W_,W_);
        W{i}{1} = W_ ; clear  W_ ;
    else
        load(param_preproc.G_filename(param_preproc.subcube,param_preproc.ch(i)),'Gw');
        for j =1:numel(Gw)
            W{i}{j,1} = any(abs(Gw{j}), 1).';
            G{i}{j,1} = Gw{j}(:,W{i}{j,1});
            Gw{j} =[];
        end
        
        % get active Fourier modes
        Wcurr = W{i}{1};
        for j = 2:numel(W{i})
            Wcurr = ((Wcurr+ W{i}{j,1})>0);
        end
        % build H and gridded data
        H{i}  = sparse(nnz(Wcurr),nnz(Wcurr));
        for j =1:numel(G{i})
            Gcurr = sparse(prod(param_nufft.No),size(G{i}{j},1));
            Gcurr(W{i}{j,1},:)= G{i}{j}'; G{i}{j} =[]; W{i}{j,1}=[];
            % gridded data
            if j==1, data{i}{1} = Gcurr*y{i}{j};
            else, data{i}{1} =  data{i}{1}  + Gcurr*y{i}{j};
            end
            Gcurr = Gcurr(Wcurr,:);
            %H
            H{i}{1} = H{i}{1} + tril(Gcurr*Gcurr'); % keep lower half only of H
        end, clear Gcurr G ;
        W{i}{1} =Wcurr; clear Wcurr;
    end
    
    data{i}{1} = data{i}{1}( W{i}{1} );


    
    fprintf('\nUpdating diag elements of H\n')
    diagh  = diag(H{i}{1});
    nnzDiagH = find(diagh);
    for idiag =1: numel(nnzDiagH)
        H{i}{1}(nnzDiagH(idiag),nnzDiagH(idiag)) = 0.5*  H{i}{1}(nnzDiagH(idiag),nnzDiagH(idiag)); 
    end, clear nnzDiagH  diagh;
    
    H{i}{1} = H{i}{1} ;%+ H{i}{1}';
   
    aW{i}{1} = 1;% no precond ...

    Sigma{i}{1} =1; % no sigma 
    
    whos H
    
    
    
    
end

end
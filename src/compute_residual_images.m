function residual_image = compute_residual_images(x, y, A, At, ...
    G, W, flag_visibility_gridding, Sigma)
% Compute residual images.
%
% Compute the residual image for each spectral channel.
%
% Parameters
% ----------
% x : double[:, :, :]
%     Wideband image cube ``[N(1), N(2), L]``.
% y : cell
%     Blocked visibilities ``{L}{nblocks}``.
% A : anonymous function
%     Measurement operator ``@[1]``.
% At : anonymous function
%     Adjoint measurement operator @[1].
% G : cell
%     Blocked gridding matrix ``{L}{nblocks}``.
% W : cell
%     Masking operators (selection of data blocks) ``{L}{nblocks}``.
% flag_visibility_gridding : bool
%     Flag to activate data dimensionality reduction via visibility gridding.
% Sigma : cell
%     Dimensionality reduction weights ``{L}{nblocks}``.
%
% Returns
% -------
% residual_image : double[:, :, :]
%     Residual image cube ``[N(1), N(2), L]``, normalised by the peak of
%     the PSF following convention in RI.
%
% -------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
% -------------------------------------------------------------------------%
%%

n_channels = size(x, 3);
residual_image = zeros(size(x));
dirac = sparse(size(x, 1)*0.5 +1, size(x, 2)*0.5 +1 ,1,size(x, 1),size(x, 2));

if flag_visibility_gridding % H  = G' +G;  G is  a lower tril matrix
    for i = 1:n_channels
        %get residual image
        Fx = A(x(:, :, i));
        r = zeros(numel(Fx), 1);
        for j = 1:length(G{i})
            if istril(G{i}{j})
                FxSlice = Fx(W{i}{j});
                res_f = Sigma{i}{j} .* (y{i}{j} - ...
                    (Sigma{i}{j} .* (G{i}{j} * FxSlice + (FxSlice' * G{i}{j})')));
                FxSlice = [];
                r(W{i}{j}) = r(W{i}{j}) + (res_f' *  G{i}{j})' + G{i}{j} *  res_f;
            else
                r(W{i}{j}) = r(W{i}{j}) + G{i}{j}' * (Sigma{i}{j} .*  ...
                    (y{i}{j} - (Sigma{i}{j} .* (G{i}{j} * Fx(W{i}{j})))));
            end
        end,  clear Fx  res_f;
        residual_image(:, :, i) = real(At(r));
        clear r;
        
        % get psf peak
        Fxdirac = A(full(dirac));
        r = zeros(numel(Fxdirac), 1);
        for j = 1:length(G{i})
            if istril(G{i}{j})
                FxSlice = Fxdirac(W{i}{j});
                res_f = (Sigma{i}{j}.^2) .* (G{i}{j} * FxSlice + (FxSlice' * G{i}{j})');
                FxSlice = [];
                r(W{i}{j}) = r(W{i}{j}) +  (res_f' *  G{i}{j})'  +  G{i}{j} *  res_f;
            else
                r(W{i}{j}) = r(W{i}{j}) + G{i}{j}' * ...
                    ((Sigma{i}{j}.^2) .* (G{i}{j} * Fxdirac(W{i}{j})));
            end
        end,  clear  Fxdirac  res_f;
        % normalise  residual image
        residual_image(:, :, i) =residual_image(:, :, i)./(max(max(real(At(r)))));
        clear r;
        
    end
else
    for i = 1:n_channels
        % get residual image
        Fx = A(x(:, :, i));
        r = zeros(numel(Fx), 1);
        for j = 1:length(G{i})
            r(W{i}{j}) = r(W{i}{j}) + G{i}{j}' * (y{i}{j} - G{i}{j} * Fx(W{i}{j}));
        end; clear Fx u2;
        residual_image(:, :, i) = real(At(r));
        clear r;
        
        % get peak of psf
        Fxdirac = A(full(dirac));
        r = zeros(numel(Fxdirac), 1);
        for j = 1:length(G{i})
            r(W{i}{j}) = r(W{i}{j}) + G{i}{j}' * (G{i}{j} * Fxdirac(W{i}{j}));
        end; clear Fxdirac u2;
        % normalise residual image
        residual_image(:, :, i) = residual_image(:, :, i)./(max(max(real(At(r)))));
        clear r;
    end
end

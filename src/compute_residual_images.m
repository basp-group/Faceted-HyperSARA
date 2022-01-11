function residual_image = compute_residual_images(x, y, A, At, ...
                                                G, W, flagDR, Sigma)
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
% flagDR : bool
%     Flag to activate DR functionality.
% Sigma : cell
%     Dimensionality reduction weights ``{L}{nblocks}``.
%
% Returns
% -------
% residual_image : double[:, :, :]
%     Residual image cube ``[N(1), N(2), L]``.
%

% -------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
% -------------------------------------------------------------------------%
%%

n_channels = size(x, 3);
residual_image = zeros(size(x));

if flagDR % H  =G' +G; and H' = H ; G is  a lower tril matrix
    for i = 1:n_channels
        Fx = A(x(:, :, i));
        r = zeros(numel(Fx), 1);
        for j = 1:length(G{i})
            if istril(G{i}{j})
             FxSlice = Fx(W{i}{j});
             res_f = Sigma{i}{j} .* (y{i}{j} - (Sigma{i}{j} .* (G{i}{j} * FxSlice + (FxSlice' * G{i}{j})'))); FxSlice = [];
                 u2 = (res_f' *  G{i}{j})'  +  G{i}{j} *  res_f; res_f = [];
        else
         % res_f = y{i}{j} - (Sigma{i}{j} .* (G{i}{j} * Fx(W{i}{j})));
             u2 = G{i}{j}' * (Sigma{i}{j} .*  (y{i}{j} - (Sigma{i}{j} .* (G{i}{j} * Fx(W{i}{j})))));
            end
        r(W{i}{j}) = r(W{i}{j}) + u2; u2 = [];
        end
        clear   Fx u2;
        residual_image(:, :, i) = real(At(r));
        clear   r;
    end
else
    for i = 1:n_channels
        Fx = A(x(:, :, i));
        r = zeros(numel(Fx), 1);
        for j = 1:length(G{i})
            % res_f = y{i}{j} - G{i}{j} * Fx(W{i}{j});
            u2 = G{i}{j}' * (y{i}{j} - G{i}{j} * Fx(W{i}{j}));
            r(W{i}{j}) = r(W{i}{j}) + u2; u2 = [];
        end; clear Fx u2;
        residual_image(:, :, i) = real(At(r));
    clear r;
    end
end

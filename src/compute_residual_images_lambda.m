function residual_image = compute_residual_images_lambda(x, y, A, At, W, ...
    apply_G, apply_Gdag)
    % Compute residual images.
    %
    % Compute the residual image for each spectral channel.
    %
    % Parameters
    % ----------
    % x : array (3d)
    %     Wideband image cube [N(1), N(2), L].
    % y : cell
    %     Blocked visibilities {L}{nblocks}.
    % G : cell
    %     Blocked gridding matrix {L}{nblocks}.
    % A : lambda
    %     Measurement operator @[1].[description]
    % At : lambda
    %     Adjoint measurement operator @[1].[description]
    % W : cell
    %     Masking operators (selection of data blocks) {L}{nblocks}.[description]
    %
    % Returns
    % -------
    % residual_image : array (3d)
    %     Residual image cube [N(1), N(2), L].
    %
    
    %-------------------------------------------------------------------------%
    %%
    % Code: P.-A. Thouvenin.
    % Last revised: [08/08/2019]
    %-------------------------------------------------------------------------%
    %%
    
    nChannels = size(x, 3);
    residual_image = zeros(size(x));
    
    for i = 1 : nChannels
        Fx = A(x(:,:,i));
        r = zeros(numel(Fx),1);
        for j = 1 : length(G{i})
            res_f = y{i}{j} - apply_G(Fx(W{i}{j}), G{i}{j});
            u2 = apply_Gdag(res_f, G{i}{j});
            r(W{i}{j}) = r(W{i}{j}) + u2;
        end
        residual_image(:,:,i) = real(At(r));
    end
    
    end
    
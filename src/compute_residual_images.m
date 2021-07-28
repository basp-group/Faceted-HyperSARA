function residual_image = compute_residual_images(x, y, G, A, At, W)
% Compute residual images.
%
% Compute the residual image for each spectral channel.
%
% Args:
%     x (array_like): wideband image cube [N(1), N(2), L].
%     y (cell): blocked visibilities {L}{nblocks}.
%     G (cell): blocked gridding matrix {L}{nblocks}.
%     A (lambda): measurement operator @[1].
%     At (lambda): adjoint measurement operator @[1].
%     W (cel): masking operators (selection of data blocks) {L}{nblocks}.
%
% Returns:
%     residual_image (array_like): residual image cube [N(1), N(2), L]
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
        res_f = y{i}{j} - G{i}{j} * Fx(W{i}{j});
        u2 = G{i}{j}' * res_f;
        r(W{i}{j}) = r(W{i}{j}) + u2;
    end
    residual_image(:,:,i) = real(At(r));
end

end

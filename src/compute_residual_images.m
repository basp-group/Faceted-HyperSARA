function residual_image = compute_residual_images(x, y, G, A, At, W)
%compute_residual_images: compute the residual images for each channel of 
% interest.
%-------------------------------------------------------------------------%
%%
% Input:
%
% > x             wideband image cube [N(1), N(2), L]
% > y             visibilities (blocked) {L}{nblocks}
% > G             gridding matrix {L}{nblocks}
% > A             measurement operator @[1]
% > At            adjoint measurement operator @[1]
% > W             masking operators (selection of data blocks) {L}{nblocks}
%
% Output:
%
% < residual_image  residual image cube [N(1), N(2), L]
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

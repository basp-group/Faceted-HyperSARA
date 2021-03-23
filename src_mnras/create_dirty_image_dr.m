function dirty_image = create_dirty_image_dr(y, At, H, T, W, Nx, Ny, No)

%! to be documented
%! create "dirty" noise matrix

nChannels = numel(y);
dirty_image = zeros(Ny, Nx, nChannels);

for l = 1:nChannels
    d_l = zeros(No, 1);
    for b = 1:numel(y{l})
        d_l(W{l}{b}) = d_l(W{l}{b}) + H{l}{b}' * (T{l}{b} .* y{l}{b});
    end
    dirty_image(:,:,l) = reshape(At(d_l),[Ny, Nx, 1]);     
end

end
    
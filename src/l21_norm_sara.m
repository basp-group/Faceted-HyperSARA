function f = l21_norm_sara(x, Psit, s)
%l21_norm_sara: compute the l21-norm of the sparsifying dictionary 
% \Psi^\dag applied to each channel of the wideband image.
%-------------------------------------------------------------------------%
%%
% Input:
%
% > x     wideband image cube [N(1), N(2), L]
% > Psit  SARA dictionary @[1]
% > s     size of the wavelet transforms (related to )
%
% Output:
%
% < f     l21-norm [1]
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% [../../2019]
%-------------------------------------------------------------------------%
%%

r = zeros(s, size(x, 3));

for l = 1:size(x, 3)
   r(:, l) = Psit(x(:, :, l)); 
end

f = sum(sqrt(sum(abs(r).^2,2)));

end

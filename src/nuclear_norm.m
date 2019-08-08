function f = nuclear_norm(x)
%nuclear_norm: compute the nuclear norm of the wideband image.
%-------------------------------------------------------------------------%
%%
% Input:
%
% > x  wideband image cube [N(1), N(2), L]
%
% Output:
%
% < f  value of the nuclear norm [1]
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
%-------------------------------------------------------------------------%
%%

[~,D,~] = svd(reshape(x, [numel(x)/size(x, 3), size(x, 3)]),'econ');
f = sum(abs(diag(D)));

end

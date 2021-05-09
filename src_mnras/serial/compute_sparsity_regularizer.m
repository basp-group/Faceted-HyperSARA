function gamma1 = compute_sparsity_regularizer(x, Psit, s, sig, alph)
% TODO: to be updated

%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [...]
%-------------------------------------------------------------------------%

r = zeros(s, size(x, 3));

for l = 1:size(x, 3)
   r(:, l) = Psit(x(:, :, l)); 
end

gamma1 = alph / (sig * sum(log(sqrt(sum(abs(r).^2,2))/sig + 1)));

end
    
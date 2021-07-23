function gamma1 = compute_sparsity_regularizer(x, Psit, s, sig, alph, regtype)
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

switch regtype
   case "log"
      gamma1 = alph / (sig * sum(log(sqrt(sum(abs(r).^2,2))/sig + 1)));
   case "inv"
      gamma1 = alph / sum(sqrt(sum(abs(r).^2,2)));
   otherwise
      gamma1 = alph / sum(sqrt(sum(abs(r).^2,2)));
end

end
    
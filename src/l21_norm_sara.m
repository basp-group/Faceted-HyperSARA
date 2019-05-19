function f = l21_norm_sara(x, Psit, s)

r = zeros(s, size(x, 3));

for l = 1:size(x, 3)
   r(:, l) = Psit(x(:, :, l)); 
end

f = sum(sqrt(sum(abs(r).^2,2)));

end

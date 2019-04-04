function x = Split_adjoint_operator_new(x_split,chunks,n,m,c)

%
x = zeros(n,m,c);
d = length(chunks);

for i = 1 : d
    x(:,:,chunks{i}) = x_split{i};
end

end

function x = Split_adjoint_operator(x_split,chunks,n,m,c)

%
x = zeros(n,m,c);

for i = 1 : length(chunks)
    x(:,:,chunks(i,1):chunks(i,2)) = x(:,:,chunks(i,1):chunks(i,2)) + x_split{i};
end

end
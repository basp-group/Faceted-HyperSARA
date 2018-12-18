function x = Split_adjoint_operator2(x_split,chunks,n,m,c)

%
x = zeros(n,m,c);

for i = 1 : length(chunks)
    x(:,:,chunks(i,1):chunks(i,2)) = x(:,:,chunks(i,1):chunks(i,2)) + x_split{i};
end

x(:,:,chunks(2,1):chunks(end-1,2)) = x(:,:,chunks(2,1):chunks(end-1,2)) / 2; 

end
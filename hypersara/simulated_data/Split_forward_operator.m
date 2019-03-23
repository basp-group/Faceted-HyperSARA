function x_split = Split_forward_operator(x,chunks)

%
x_split = cell(length(chunks),1);

for i = 1 : length(chunks)
    x_split{i} = x(:,:,chunks(i,1):chunks(i,2));
end

end
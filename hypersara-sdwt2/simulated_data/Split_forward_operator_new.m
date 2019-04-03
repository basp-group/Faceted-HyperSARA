function x_split = Split_forward_operator_new(x,chunks)

%
d = length(chunks);
x_split = cell(d,1);

for i = 1 : d
    x_split{i} = x(:,:,chunks{i});
end

end

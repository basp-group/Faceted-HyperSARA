function x_split = Split_forward_operator(x,I,dims,Q)

%
x_split = cell(Q,1);

for q = 1:Q
    x_split{q} = x(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2),:);
end

end
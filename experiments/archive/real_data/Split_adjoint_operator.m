function x = Split_adjoint_operator(x_split,I,dims,Q,n,m,c)

%
x = zeros(n,m,c);

for q = 1:Q
    x(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2),:) = x_split{q};
end

end



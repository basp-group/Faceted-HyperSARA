function v = Split_sparsity_forward_operator(x,Sp,Psit)

%
x_split = Sp(x);

for  i = 1 : length(x_split)
   v{i} = Psit(x_split{i});
end

end
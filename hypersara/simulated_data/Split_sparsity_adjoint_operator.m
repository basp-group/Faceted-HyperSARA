function x = Split_sparsity_adjoint_operator(v,Spt,Psi)

%
for i = 1 : length(v)
x_split{i} = Psi(v{i});
end

x = Spt(x_split);

end
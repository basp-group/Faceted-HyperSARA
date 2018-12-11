function x = Croping_adjoint_operator(s,Mt,c,Ny,Nx)

for i = 1 : c
    x(:,:,i) = reshape(Mt * s{i},Ny,Nx);
end

end
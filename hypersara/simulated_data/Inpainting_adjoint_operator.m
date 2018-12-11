function x = Inpainting_adjoint_operator(s,Lt_inp,flag_ind,c,Ny,Nx)

empt = zeros(size(Lt_inp{1}{1} * s{flag_ind{1}(1)}));

for i = 1 : c
    temp = empt;
    for j = 1 : length(Lt_inp{i})
        temp = temp + Lt_inp{i}{j} * s{flag_ind{i}(j)};
        x(:,:,i) = reshape(temp, Ny,Nx);
    end
end

end
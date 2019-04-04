function x = Inpainting_adjoint_operator_crop(s,Lt_inp,flag_ind,c)

empt = zeros(size(Lt_inp{1}{1} * s{flag_ind{1}(1)}));

for i = 1 : c
    x{i} = empt;
    for j = 1 : length(Lt_inp{i})
        x{i} = x{i} + Lt_inp{i}{j} * s{flag_ind{i}(j)};
    end
end

end
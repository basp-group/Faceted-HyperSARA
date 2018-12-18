function s = Inpainting_forward_operator_crop(x,L_inp,flag_ind,c)

empt = zeros(size(L_inp{1}{1} * x{1}));

for i = 1 : c
    s{i} = empt;
    for j = 1 : length(L_inp{i})
        s{i} = s{i} + L_inp{i}{j} * x{flag_ind{i}(j)};
    end
end

end
function s = Inpainting_forward_operator(x,L_inp,flag_ind,c)

temp = x(:,:,1);
empt = zeros(size(L_inp{1}{1} * temp(:)));

for i = 1 : c
    s{i} = empt;
    for j = 1 : length(L_inp{i})
        temp = x(:,:,flag_ind{i}(j));
        s{i} = s{i} + L_inp{i}{j} * temp(:);
    end
end

end
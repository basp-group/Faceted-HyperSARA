function s = Croping_forward_operator(x,M,c)

for i = 1 : c
    temp = x(:,:,i);
    s{i} = M * temp(:);
end

end
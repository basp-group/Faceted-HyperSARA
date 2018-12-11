function s = Croping_forward_operator_crop(x,M,c)

for i = 1 : c
    s{i} = M * x{i};
end

end
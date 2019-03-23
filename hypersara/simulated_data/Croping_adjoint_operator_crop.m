function x = Croping_adjoint_operator_crop(s,Mt,c)

for i = 1 : c
    x{i} = Mt * s{i};
end

end
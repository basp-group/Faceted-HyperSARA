function im = placeFacet(PsiSty, I_overlap, dims_overlap, N)

start = 1;
im = zeros(N);
for m = 1:size(I_overlap, 1)
    corner = I_overlap(m, :);     % index corner after extension
    segSize = dims_overlap(m, :); % facet dimension after extension 
    s = prod(dims_overlap(m, :));
    im(corner(1)+1:corner(1)+segSize(1),corner(2)+1:corner(2)+segSize(2)) = ...
    reshape(PsiSty(start:start+s-1), dims_overlap(m, :));
    start = start + s;
end
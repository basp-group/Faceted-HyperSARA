function im = place2DSegment(im, PsiSty, I_facet_min, dims_facet_max)

% input min I_overlap + max dims_overlap

im(I_facet_min(1)+1:I_facet_min(1)+dims_facet_max(1),I_facet_min(2)+1:I_facet_min(2)+dims_facet_max(2), :) = ...
im(I_facet_min(1)+1:I_facet_min(1)+dims_facet_max(1),I_facet_min(2)+1:I_facet_min(2)+dims_facet_max(2), :) + PsiSty;


end
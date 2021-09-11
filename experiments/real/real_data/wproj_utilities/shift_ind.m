function ind = shift_ind(ind, Nx, Ny)

i = ind(:, 1);
j = ind(:, 2);
% size assumed even
i_ = zeros(length(i), 1);
j_ = zeros(length(j), 1);
i_(i < Nx / 2 + 1) = i(i < Nx / 2 + 1) + Nx / 2;
i_(i > Nx / 2) = i(i > Nx / 2) - Nx / 2;
j_(j < Ny / 2 + 1) = j(j < Ny / 2 + 1) + Ny / 2;
j_(j > Ny / 2) = j(j > Ny / 2) - Ny / 2;

ind = [i_(:) j_(:)];

end

function [Iq, dims_q, dims_oq, dims_overlap_ref_q, I_overlap_q, ...
    dims_overlap_q, status_q, offsetLq, offsetRq, Ncoefs_q, temLIdxs_q, ...
    temRIdxs_q, overlap_g_south, overlap_g_east, overlap_g_south_east, ...
    overlap, w, crop_low_rank, crop_sparsity] = fhs_setup_priors(Qx, Qy, I, dims, ...
    dims_o, dims_overlap_ref, I_overlap, dims_overlap, status, offsetL, ...
    offsetR, Ncoefs, temLIdxs, temRIdxs, window_type, overlap_size)
% Setup all the composite arrays for the faceted low-rankness and average
% joint sparsity priors involed in Faceted HyperSARA.
%
% Parameters
% ----------
% Qx : int
%     Number of facets (axis x).
% Qy : int
%     Number of facets (axis y).
% I : int[:]
%     Starting index of the non-overlapping tile ``[1, 2]``.
% dims : int[:]
%     Size of the non-overlapping tile ``[1, 2]``.
% dims_o : int[:]
%     Dimension of a facet (with overlap) ``[1, 2]``.
% dims_overlap_ref : int[:]
%     Dimension of the facet ``[1, 2]``.
% I_overlap : int[:]
%     Starting index of the overlapping facets ``[Q, 1]``.
% dims_overlap : int[:, :]
%     Size of the overlapping facets ``[Q, 2]``.
% status : int[:]
%     Status of the current facet (last or first facet along vert. / hrz.
%     direction) ``[1, 2]``.
% offsetL : int[:, :]
%     Amount of zero-pading from the "left" ``[Q, 2]``.
% offsetR : int[:, :]
%     Amount of zero-padding from the "right" ``[Q, 2]``.
% Ncoefs : cell of int[:, :]
%     Size of the wavelet decomposition across each scale.
% temLIdxs : int[:, :]
%     Amount of cropping from the "left" ``[Q, 2]``.
% temRIdxs : int[:, :]
%     Amount of cropping from the "right" ``[Q, 2]``.
% window_type : string
%     Type of apodization window considered for the faceted low-rankness
%     prior.
% overlap_size : int[:]
%     Number of overlapping pixels in the current facet along each
%     direction ``[1, 2]``.
%
% Returns
% -------
% Iq : int[:]
%     Starting index of the non-overlapping base facet ``[1, 2]``.
% dims_q : int[:]
%     Dimensions of the non-overlapping base facet ``[1, 2]``.
% dims_oq : int[:]
%     Diemnsions of the overlapping facet ``[1, 2]``.
% dims_overlap_ref_q : int[:]
%     Dimension of the facet ``[1, 2]``.
% I_overlap_q : int[:]
%     Starting index of the facet ``[1, 2]``.
% dims_overlap_q : int[:]
%     Dimensions of the facet ``[1, 2]``.
% status_q : int[:]
%     Status of the current facet (last or first facet along vert. or hrz.
%     direction) [ndict, 2].
% offsetLq : int[:]
%     Amount of zero-pading from the "left" ``[1, 2]``.
% offsetRq : int[:]
%     Amount of zero-pading from the "right" ``[1, 2]``.
% Ncoefs_q : int[:, :]
%     Size of the wavelet decompositions at each scale.
% temLIdxs_q : int[:]
%     Amount of cropping from the "left" ``[1, 2]``.
% temRIdxs_q : int[:]
%     Amount of cropping from the "right" ``[1, 2]``.
% overlap_g_south : int[:]
%     Size of the overlap for the south neighbour [2,1].
% overlap_g_east : int[:]
%     Size of the overlap for the east neighbour [2,1].
% overlap_g_south_east : int[:]
%     Size of the overlap for the south-east neighbour [2,1].
% overlap : int[:]
%     Number of overlapping pixels along each direction.
% w : double[:, :]
%     Apodization window considered for the faceted low-rankness prior.
% crop_low_rank : int[:]
%     [Relative cropping necessary for the faceted low-rankness prior
%     ``[1, 2]``.
% crop_sparsity : int[:]
%     Relative cropping necessary for the faceted joint-sparsity prior
%     ``[1, 2]``.
%

% define composite variables (local to a given worker)
% /!\ only simple indexing allowed into Composite objects from the master
% node
Iq = Composite();
dims_q = Composite();
temLIdxs_q = Composite();
temRIdxs_q = Composite();
I_overlap_q = Composite();
dims_overlap_q = Composite();
dims_overlap_ref_q = Composite();
status_q = Composite();
Ncoefs_q = Composite();
offsetLq = Composite();
offsetRq = Composite();
% dimension of the ghost cells
overlap_g_south = Composite();
overlap_g_east = Composite();
overlap_g_south_east = Composite();
overlap = Composite();
% constant overlap: facet size
dims_oq = Composite();
crop_low_rank = Composite();
crop_sparsity = Composite();

Q = Qx * Qy;

flag_overlap = any(overlap_size > 0);

% initialize composite variables and constants for the faceted prior (l21
% and nuclear norms)
for q = 1:Q
    Iq{q} = I(q, :);
    dims_q{q} = dims(q, :);
    temLIdxs_q{q} = temLIdxs{q};
    temRIdxs_q{q} = temRIdxs{q};
    I_overlap_q{q} = I_overlap{q};
    dims_overlap_q{q} = dims_overlap{q};
    status_q{q} = status(q, :);
    Ncoefs_q{q} = Ncoefs{q};

    % additional composite variables (for the zero padding / boundary conditions)
    dims_overlap_ref_q{q} = dims_overlap_ref(q, :);
    offsetLq{q} = offsetL(q, :);
    offsetRq{q} = offsetR(q, :);

    % check which overlap is the largest (sdwt2 or nuclear norms)
    bool_crop = dims_overlap_ref(q, :) >= dims_o(q, :); % 0: nuclear norm largest, 1: dwt2 largest
    tmp_crop_l21 = [0, 0];
    tmp_crop_nuclear = [0, 0];
    if bool_crop(1)
        % sdwt2 largest
        tmp_crop_nuclear(1) = dims_overlap_ref(q, 1) - dims_o(q, 1);
    else
        tmp_crop_l21(1) = dims_o(q, 1) - dims_overlap_ref(q, 1);
    end

    if bool_crop(2)
        % sdwt2 largest
        tmp_crop_nuclear(2) = dims_overlap_ref(q, 2) - dims_o(q, 2);
    else
        tmp_crop_l21(2) = dims_o(q, 2) - dims_overlap_ref(q, 2);
    end
    crop_sparsity{q} = tmp_crop_l21;
    crop_low_rank{q} = tmp_crop_nuclear;
    overlap{q} = max(max(dims_overlap{q}) - dims(q, :), dims_o(q, :) - dims(q, :));
    % amount of overlap necessary for each facet
    dims_oq{q} = dims_o(q, :);
end

% define apodization window for the faceted nuclear norm prior
w = Composite();
for q = 1:Q
    [qy, qx] = ind2sub([Qy, Qx], q);
    if qy < Qy
        % S (qy+1, qx)
        overlap_g_south{q} = overlap{(qx - 1) * Qy + qy + 1};

        if qx < Qx
            % SE (qy+1, qx+1)
            overlap_g_south_east{q} = overlap{qx * Qy + qy + 1};
        else
            overlap_g_south_east{q} = [0, 0];
        end
    else
        overlap_g_south{q} = [0, 0];
        overlap_g_south_east{q} = [0, 0];
    end
    if qx < Qx
        % E (qy, qx+1)
        overlap_g_east{q} = overlap{qx * Qy + qy};
    else
        overlap_g_east{q} = [0, 0];
    end

    % define the weights (depends on the position of the facet inside the grid)
    if flag_overlap
        w{q} = generate_weights(qx, qy, Qx, Qy, window_type, dims(q, :), dims_o(q, :), overlap_size);
    else
        w{q} = ones(dims_o(q, :)); % ! in this case, dims = dims_o
    end
end

end

function [Iq, dims_q, dims_oq, dims_overlap_ref_q, I_overlap_q, ...
    dims_overlap_q, status_q, offsetLq, offsetRq, Ncoefs_q, temLIdxs_q, ...
    temRIdxs_q, overlap_g_south, overlap_g_east, overlap_g_south_east, overlap, ...
    w, crop_nuclear, crop_l21] = setup_priors(Qx, Qy, I, dims, dims_o, ...
    dims_overlap_ref, I_overlap, dims_overlap, status, offsetL, offsetR, ...
    Ncoefs, temLIdxs, temRIdxs, window_type, d)

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
crop_nuclear = Composite();
crop_l21 = Composite();

Q = Qx*Qy;

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
    dims_overlap_ref_q{q} = dims_overlap_ref(q,:);
    offsetLq{q} = offsetL(q,:);
    offsetRq{q} = offsetR(q,:);
    
    % check which overlap is the largest (sdwt2 or nuclear norms)
    bool_crop = dims_overlap_ref(q,:) >= dims_o(q,:); % 0: nuclear norm largest, 1: dwt2 largest
    tmp_crop_l21 = [0, 0];
    tmp_crop_nuclear = [0, 0];
    if bool_crop(1)
        % sdwt2 largest
        tmp_crop_nuclear(1) = dims_overlap_ref(q,1) - dims_o(q,1);
    else
        tmp_crop_l21(1) = dims_o(q,1) - dims_overlap_ref(q,1);
    end
    
    if bool_crop(2)
        % sdwt2 largest
        tmp_crop_nuclear(2) = dims_overlap_ref(q,2) - dims_o(q,2);
    else
        tmp_crop_l21(2) = dims_o(q,2) - dims_overlap_ref(q,2);
    end
    crop_l21{q} = tmp_crop_l21;
    crop_nuclear{q} = tmp_crop_nuclear;    
    overlap{q} = max(max(dims_overlap{q}) - dims(q,:), dims_o(q, :) - dims(q,:)); 
    % amount of overlap necessary for each facet
    dims_oq{q} = dims_o(q, :);
end

% define apodization window for the faceted nuclear norm prior
w = Composite();
for q = 1:Q
    [qy, qx] = ind2sub([Qy, Qx], q);
    if qy < Qy
        % S (qy+1, qx)
        overlap_g_south{q} = overlap{(qx-1)*Qy + qy+1};
        
        if qx < Qx
            % SE (qy+1, qx+1)
            overlap_g_south_east{q} = overlap{qx*Qy + qy+1};
        else
            overlap_g_south_east{q} = [0, 0];
        end
    else
        overlap_g_south{q} = [0, 0];
        overlap_g_south_east{q} = [0, 0];
    end
    if qx < Qx
        % E (qy, qx+1)
        overlap_g_east{q} = overlap{qx*Qy + qy};
    else
        overlap_g_east{q} = [0, 0];
    end
    
    % define the weights (depends on the position of the facet inside the grid)
    w{q} =  generate_weights(qx, qy, Qx, Qy, window_type, dims(q,:), dims_o(q,:), [d,d]);
end

end
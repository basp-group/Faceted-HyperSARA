function val = op_norm_sdwt2(N, I, dims, nlevel, wavelet, L, tol, max_iter, verbose)
%% computes the maximum eigen value of the compund operator At*A
y = randn(N);
y = y/norm(y(:));
init_val = 1;

Q = size(dims, 1);
SPsitLx = cell(Q, 1);
PsiStu = cell(Q, 1);
Ij = cell(Q, 1);
dims_PsitLx = cell(Q, 1);
Ncoefs = cell(Q, 1);
[~, ~, I_overlap_ref, dims_overlap_ref, I_overlap, dims_overlap, ...
    I_overlap_nc, dims_overlap_nc, status, offset, ~, ~, offsetL, offsetR] = generate_segdwt_indices(N, I, dims, nlevel, wavelet, L);

for k = 1:max_iter
    % x = A(y);
    for q = 1:Q
        zerosNum = dims_overlap_ref(q,:) + offsetL(q,:) + offsetR(q,:);
        x_overlap = zeros([zerosNum(:)',1]);
        
        x_overlap(offsetL(q,1)+1:end-offsetR(q,1),...
            offsetL(q,2)+1:end-offsetR(q,2))...
            = y(I_overlap_ref(q, 1)+1:I_overlap_ref(q, 1)+dims_overlap_ref(q, 1), ...
            I_overlap_ref(q, 2)+1:I_overlap_ref(q, 2)+dims_overlap_ref(q, 2));
        
        % forward operator [put the following instructions into a parfeval for parallelisation]
        [SPsitLx{q}, Ij{q}, dims_PsitLx{q}, Ncoefs{q}] = sdwt2_sara(x_overlap, I(q, :), dims(q, :), offset, status(q,:), nlevel, wavelet);
        
        % inverse operator (for a single facet) (inverse = adjoin for spd, properly implement the adjoint operator for different boundary conditions)
        PsiStu{q} = isdwt2_sara(SPsitLx{q}, I(q, :), dims(q, :), I_overlap_nc{q}, dims_overlap_nc{q}, Ncoefs{q}, N, nlevel, wavelet);
    end

    % y = At(x);
    y = zeros(N);
    for q = 1:Q
        y = place2DSegment(y, PsiStu{q}, I_overlap{q}, dims_overlap{q});
    end

    val = norm(y(:));
    rel_var = abs(val-init_val)/init_val;
    
    if (verbose > 1)
        fprintf('Iter = %i, norm = %e \n', k, val);
    end
    if (rel_var < tol)
       break;
    end
    
    init_val = val;
    y = y/val;  
end

if (verbose > 0)
    fprintf('Norm = %e \n\n', val);
end

end


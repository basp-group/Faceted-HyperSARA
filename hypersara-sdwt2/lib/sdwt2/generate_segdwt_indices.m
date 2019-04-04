function [I_overlap_ref_nc, dims_overlap_ref_nc, I_overlap_ref, ...
    dims_overlap_ref, I_overlap, dims_overlap, I_overlap_nc, ...
    dims_overlap_nc, status, offset, rJnew, rJnew_ref, offsetL, offsetR] = generate_segdwt_indices(N, I, dims, J, wavelet, L)
%%
M = numel(wavelet);
rJnew_ref = (2^J-1)*(max(L(:))-2);
rJnew = zeros(M, 1);
for m = 1:M
    if ~strcmp(wavelet{m}, 'self')
        rJnew(m) = (2^J-1)*(L(m) - 2); % r_red(J) (4.12)
    else
        rJnew(m) = 0; % do not forget to add mod(I(q, :), 2^J); (see l.64 below) to the normal offset
    end
end
offset = rJnew_ref - rJnew;
dim = 2;
Q = size(I, 1);

I_overlap = cell(Q, 1);
I_overlap_nc = cell(Q, 1);
I_overlap_ref = zeros(Q, 2);
I_overlap_ref_nc = zeros(Q, 2);
status = zeros(Q, 2);
dims_overlap = cell(Q, 1);
dims_overlap_nc = cell(Q, 1);
dims_overlap_ref = zeros(Q, 2);
dims_overlap_ref_nc = zeros(Q, 2);
offsetL = zeros(Q, 2);
offsetR = zeros(Q, 2);

for q = 1:Q                            % define facets, read appropriate portion of the image of interest
    
    LnoRrows = rJnew_ref + mod(I(q,1),2^J);        % entension width
    LnoRcols = rJnew_ref + mod(I(q,2),2^J);        % idem
    I_overlap_ref_nc(q,1)= I(q,1)-LnoRrows;       % starting index after left extension
    I_overlap_ref_nc(q,2)= I(q,2)-LnoRcols;       % idem
    dims_overlap_ref_nc(q,1)= dims(q,1)+LnoRrows; % dimension after left extension
    dims_overlap_ref_nc(q,2)= dims(q,2)+LnoRcols; % idem
    
    % note: extension width and dimensions are amended later, depending on the
    % position of the facet (first/last along each dimension)
    
    dimensions = dims_overlap_ref_nc(q, :);
    corner = I_overlap_ref_nc(q, :);
    
    for i=1:dim
        if corner(i)<0
            offsetL(q,i)= -corner(i);
            dimensions(i) = dimensions(i) + corner(i);
            corner(i)=0;
        end
        if corner(i)+dimensions(i)>=N(i)
            offsetR(q,i) = corner(i) + dimensions(i) - N(i);
            dimensions(i) = N(i)-corner(i);
        end
    end
    I_overlap_ref(q, :) = corner;
    dims_overlap_ref(q, :) = dimensions;
    
    for i=1:dim
        bool_first = false;
        if I(q, i)==0
            status(q,i) = -1; %'first';
            offsetL(q,i) = 0;
            bool_first = true;
        end
        if I(q, i)+dims(q, i) == N(i)
            if bool_first
                status(q,i) = NaN;
            else
                status(q,i) = 1; %'last';
            end
            offsetR(q,i) = 0;
        end
    end
    
    I_overlap{q} = zeros(M, 2);
    dims_overlap{q} = zeros(M, 2);
    I_overlap_nc{q} = zeros(M, 2);
    dims_overlap_nc{q} = zeros(M, 2);
    
    for m = 1:M
        if strcmp(wavelet{m}, 'self')
            I_overlap_nc{q}(m, :) = I(q, :); % no overlap for the dirac basis
            dims_overlap_nc{q}(m, :) = dims(q, :);
            I_overlap{q}(m, :) = I(q, :);
            dims_overlap{q}(m, :) = dims(q, :);
        else
            LnoRrows = rJnew(m) + mod(I(q,1),2^J);         % extension width
            LnoRcols = rJnew(m) + mod(I(q,2),2^J);         % idem
            I_overlap_nc{q}(m, 1) = I(q,1)-LnoRrows;       % starting index after left extension
            I_overlap_nc{q}(m, 2) = I(q,2)-LnoRcols;       % idem
            dims_overlap_nc{q}(m, 1) = dims(q,1)+LnoRrows; % dimension after left extension
            dims_overlap_nc{q}(m, 2) = dims(q,2)+LnoRcols; % idem         
            
            corner = I_overlap_nc{q}(m,:);
            dimensions = dims_overlap_nc{q}(m,:);
            for i=1:dim
                if corner(i)<0
                    dimensions(i) = dimensions(i) + corner(i);
                    corner(i)=0;
                end
                if corner(i)+dimensions(i)>=N(i)
                    dimensions(i) = N(i)-corner(i);
                end
            end
            I_overlap{q}(m,:) = corner;
            dims_overlap{q}(m,:) = dimensions;
        end
    end
end

end

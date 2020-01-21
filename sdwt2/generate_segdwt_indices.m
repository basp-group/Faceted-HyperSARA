function [I_overlap_ref_nc, dims_overlap_ref_nc, I_overlap_ref, ...
    dims_overlap_ref, I_overlap, dims_overlap, I_overlap_nc, ...
    dims_overlap_nc, status, offset, offsetL, offsetR, Ncoefs, temLIdxs, temRIdxs] = generate_segdwt_indices(N, I, dims, J, wavelet, L)
%
% 1.to be purged from unused lines of codes (check with the article again)
% 2.possibly simplify the loops in this function, write consistent comments
% here (once acceleration is confirmed)
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
temLIdxs = cell(Q, 1);
temRIdxs = cell(Q, 1);

% determine if the Dirac basis ('self') is among the dictionaries
dirac_present = any(ismember(wavelet, 'self'));
if dirac_present
   s_Ncoefs = (M-1)*(J+1)+1; % number of rows in the Ncoefs matrix (sizes for all the local decomposition levels)
else
   s_Ncoefs = M*(J+1);
end
Ncoefs = cell(Q, 1); % number of coefficients corresponding to each dictionary, at each scale of interest

for q = 1:Q % define facets, read appropriate portion of the image of interest
    
    LnoRrows = rJnew_ref + mod(I(q,1),2^J);       % entension width
    LnoRcols = rJnew_ref + mod(I(q,2),2^J);       % idem
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
    
    % try to simplify this part...
    for i=1:dim
        bool_first = false;
        if I(q, i)==0
            status(q,i) = -1; % first
            offsetL(q,i) = 0;
            bool_first = true;
        end
        if I(q, i)+dims(q, i) == N(i)% add loop over the type of dictionary
            if bool_first % first & last
                status(q,i) = NaN;
            else
                status(q,i) = 1; % last
            end
            offsetR(q,i) = 0;
        end
    end
    
    % Compute starting index/size of the overlapping facets
    I_overlap{q} = zeros(M, 2);
    dims_overlap{q} = zeros(M, 2);
    I_overlap_nc{q} = zeros(M, 2);
    dims_overlap_nc{q} = zeros(M, 2);
    temLIdxs{q} = zeros(M, 2); % offset for left cropping in isdwt2 (see possible simplification)
    temRIdxs{q} = zeros(M, 2); % offset for right cropping in isdwt2 (see possible simplification)
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
                    temLIdxs{q}(m,i)= -corner(i);
                    dimensions(i) = dimensions(i) + corner(i);
                    corner(i)=0;
                end
                if corner(i)+dimensions(i)>=N(i)
                    temRIdxs{q}(m,i) = corner(i)+dimensions(i) - N(i);
                    dimensions(i) = N(i)-corner(i);
                end
            end
            I_overlap{q}(m,:) = corner;
            dims_overlap{q}(m,:) = dimensions;
        end
    end
    
    % Compute number of coefficients for each dictionary at each scale
    Ncoefs{q} = zeros(s_Ncoefs,2);
    Sn = I(q, :);
    Snplus1 = dims(q,:) + Sn;
    id_Ncoefs = 0;
    for m = 1:M
        if ~strcmp(wavelet{m}, 'self')
            for d = 1:dim
                Snj = floor(Sn(d)./(2.^(1:J).'));
                if status(q, d) > 0 || isnan(status(q, d)) % last / first & last
                    Snplus1j = floor(2.^(-(1:J).').*Snplus1(d)+(1-2.^(-(1:J).'))*(L(m)-1));
                else
                    Snplus1j = floor(Snplus1(d)./(2.^(1:J).'));
                end
                Ncoefs{q}(id_Ncoefs+1:id_Ncoefs+J,d) = Snplus1j - Snj; % [P.-A.] (4.19) Ncoefs{q}((m-1)*(J+1)+1:m*(J+1)-1,d)
            end
%             Ncoefs{q}(m*(J+1),:) = Ncoefs{q}(m*(J+1)-1,:); % change size of Ncoefs for the self dictionary (just write dims..., no need for any overlap)
            Ncoefs{q}(id_Ncoefs+J+1,:) = Ncoefs{q}(id_Ncoefs+J,:);
            id_Ncoefs = id_Ncoefs + (J+1);
        else
            Ncoefs{q}(id_Ncoefs+1, :) = dims(q, :);
            id_Ncoefs = id_Ncoefs + 1;
        end
    end
end

end

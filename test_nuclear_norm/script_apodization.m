clc; clear all; close all;
format compact;

addpath ../sdwt2

%% Check apodization functions, see their use with overlapping facets 
% (in the reconstruction process)

N = [512, 512];
Qx = 3;
Qy = 3;
Q = Qx*Qy;
tol = 1e-3;

% regular image facets (no overlap)
rg_y = domain_decomposition(Qy, N(1));
rg_x = domain_decomposition(Qx, N(2));
I = zeros(Q, 2);
dims = zeros(Q, 2);
for qx = 1:Qx
    for qy = 1:Qy
        q = (qx-1)*Qy+qy;
        I(q, :) = [rg_y(qy, 1)-1, rg_x(qx, 1)-1];
        dims(q, :) = [rg_y(qy,2)-rg_y(qy,1)+1, rg_x(qx,2)-rg_x(qx,1)+1];
    end
end

% w = ones(N);
% % N = Npadded - Nfacet + 1;
% 
% % Generate the 2D Gaussian window
% sigx = 5e-1;
% sigy = 5e-1;
% alpha_x = 2*(N(2) - 1)*sigx;
% alpha_y = 2*(N(1) - 1)*sigy;
% wy = gausswin(N(1), alpha_y);
% wx = gausswin(N(2), alpha_x); % possibly take a different window...

% wx = kaiser(N(2));
% wx(N/2+1:end) = 1;
% wy = kaiser(N(1));
% wy(N/2+1:end) = 1;
% wg = wy.*(wx.');

% wx2 = kaiser(N(2));
% wx2(1:N/2) = 1;
% wy2 = kaiser(N(1));
% wy2(1:N/2) = 1;
% wg2 = wy2.*(wx2.');
% 
% wsum = wg + wg2;
% wg = wg./wsum;
% wg2 = wg2./wsum;

% % Convolve with a unit window of the facet size
% wg = wy.*(wx.');
% w = ones(N);
% w = conv2(wg, w, 'same');
% 
% % Threshold the smallest coefficients to zero / renormalize the window
% w(w < 1e-5) = 1e-5; % set threshold as a parameter
% w = w/max(w(:)); % normalized window

%% Option 2
% just linear interpolation (needs to be checked)

% maximum overlap in sdwt2, d can be taken smaller in practice
d = 50; % (2^nlevel - 1)*(filter_size-1);
rg_yo = domain_decomposition_overlap2(Qy, N(1), d);
rg_xo = domain_decomposition_overlap2(Qx, N(2), d);
Io = zeros(Q, 2);
dims_o = zeros(Q, 2);
for qx = 1:Qx
    for qy = 1:Qy
        q = (qx-1)*Qy+qy;
        Io(q, :) = [rg_yo(qy, 1)-1, rg_xo(qx, 1)-1];
        dims_o(q, :) = [rg_yo(qy,2)-rg_yo(qy,1)+1, rg_xo(qx,2)-rg_xo(qx,1)+1];
    end
end

% create weight matrix Wo
Wo = zeros(N);
for q = 1:Q
   Wo(Io(q,1)+1:Io(q,1)+dims_o(q,1), Io(q,2)+1:Io(q,2)+dims_o(q,2)) = ...
       Wo(Io(q,1)+1:Io(q,1)+dims_o(q,1), Io(q,2)+1:Io(q,2)+dims_o(q,2)) + ones(dims_o(q,:)); 
end
Wo = 1./Wo;

% 
wd = [linspace(tol, 1-tol, d), ones(1, N(1)-2*d), linspace(1-tol, tol, d)];
wd2 = wd.*(wd'); % not bad

% test sum 
% hrz
test_h = wd2(:,end-d+1:end) + wd2(:,1:d);
% vert
test_v = wd2(end-d+1:end,:) + wd2(1:d,:);
% diag
test_d = wd2(1:d, 1:d) + wd2(end-d+1:end, end-d+1:end); % ok

%% option 2
% V = zeros(6, 6);
% V(:,1) = tol;
% V(:,[2,3]) = [tol; 1-tol; 1; 1; 1-tol; tol].*ones(1, 2);
% V(:, 4:6) = fliplr(V(:, 1:3));
% 
% x = [1, d, d+1, N(2)-d, N(2)-d+1, N(2)];
% y = [1, d, d+1, N(1)-d, N(1)-d+1, N(1)];
% [X, Y] = meshgrid(x, y);
% [Xq, Yq] = meshgrid(1:N(2), 1:N(1));
% Vq = interp2(X, Y, V, Xq, Yq, 'linear');
% 
% % true bilinear interpolation (to be kept in the final algo?)
% test = Vq(1:d,1:d) + Vq(end-d+1:end, end-d+1:end) + ...
%        Vq(1:d,end-d+1:end) + Vq(end-d+1:end, 1:d); 
% 
% % test with multiple facets, constant overlap
% w = cell(Q, 1);
% for q = 1:Q
%    [qy, qx] = ind2sub([Qy, Qx], q);
%    Vq = V;
%     if qx == 1
%         Vq(1:end-1, 1:2) = 1;
%         if qy == 1
%             Vq(1:2,1:end-1) = 1;
%         elseif qy == Qy
%             Vq(end-1:end, :) = 1;
%         end
%     elseif qx == Qx
%         Vq(1:end-2, end-1:end) = 1;
%         if qy == 1
%             Vq(1:2,3:end) = 1;
%         elseif qy == Qy
%             Vq(end-1:end, :) = 1;
%             
%         end
%     end
%     xx = [1, d, d+1, dims_o(q,2)-d, dims_o(q,2)-d+1, dims_o(q,2)];
%     yy = [1, d, d+1, dims_o(q,1)-d, dims_o(q,1)-d+1, dims_o(q,1)];
%     [XX, YY] = meshgrid(xx, yy);
%     [Xq, Yq] = meshgrid(1:dims_o(q,2), 1:dims_o(q,1));
%     w{q} = interp2(XX, YY, Vq, Xq, Yq, 'linear');     
% end

% ok, perfect... use this procedure instead (easier to implement...)
w = cell(Q, 1);
for q = 1:Q
   [qy, qx] = ind2sub([Qy, Qx], q);
    if qx == 1
        wdx = [ones(1, d), ones(1, dims_o(q,2)-2*d), linspace(1-tol, tol, d)];
    elseif qx == Qx
        wdx = [linspace(tol, 1-tol, d), ones(1, dims_o(q,2)-2*d), ones(1, d)];
    else
        wdx = [linspace(tol, 1-tol, d), ones(1, dims_o(q,2)-2*d), linspace(1-tol, tol, d)];
    end
    
    if qy == 1
        wdy = [ones(1, d), ones(1, dims_o(q,1)-2*d), linspace(1-tol, tol, d)];
    elseif qy == Qy
        wdy = [linspace(tol, 1-tol, d), ones(1, dims_o(q,1)-2*d), ones(1, d)];
    else
        wdy = [linspace(tol, 1-tol, d), ones(1, dims_o(q,1)-2*d), linspace(1-tol, tol, d)];
    end
    w{q} = (wdy.').*wdx;    
end

figure;
for q = 1:Q
    [qy, qx] = ind2sub([Qy, Qx], q);
    qm = (qy-1)*Qx + qx;
    subplot(Qy, Qx, qm);
    imagesc(w2{q}, [0, 1]); % colormap gray; 
    axis off; axis image; colorbar % to be checked again (problem in the way the images are displayed)
end

test = w{1}(end-d+1:end, end-d+1:end) + w{2}(1:d, end-d+1:end) ...
       + w{4}(end-d+1:end, 1:d) + w{5}(1:d, 1:d);

%% Different apodization function (smooth functions, choice to be possibly optimized)

window_name = 'hamming'; 
wc = window(window_name,2*d).'; % make sure I just generate the part I'm interested in
wr = window(window_name,2*d).';

wc2 = [wc(1:d), ones(1, N(2)-d)];
wr2 = [wr(1:d), ones(1, N(2)-d)]; 

% [maskr,maskc] = meshgrid(wr2,wc2);
% w = maskr.*maskc;
% w2 = (wr2.').*wc2; % equivalent to the above two instructions (and better)

w2 = cell(Q, 1);
for q = 1:Q
   [qy, qx] = ind2sub([Qy, Qx], q);
   dims_diff = dims_o(q, :) - dims(q, :);
   if qx > 1
       wc = window(window_name,2*dims_diff(2)).'; % make sure I just generate the part I'm interested in
       wc = [wc(1:dims_diff(2)), ones(1, dims(q, 2))];
   else
       wc = ones(1, dims_o(q, 2));
   end
   if qy > 1
       wr = window(window_name,2*dims_diff(1)).';
       wr =  [wr(1:dims_diff(1)), ones(1, dims(q, 1))];
   else
       wr = ones(1, dims_o(q, 1));
   end
   w2{q} = (wr.').*wc;
end

figure;
min_w = Inf;
for q = 1:Q
    [qy, qx] = ind2sub([Qy, Qx], q);
    qm = (qy-1)*Qx + qx;
    subplot(Qy, Qx, qm);
    imagesc(w2{q}, [0, 1]); axis off; axis image; % to be checked again (problem in the way the images are displayed)
    colorbar
    min_w = min(min_w, min(w2{q}(:)));
end
% keep hamming window for the moment, change the window later on
%


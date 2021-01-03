% Illustration forthe weights
clc; clear all; close all;

addpath ../lib/utils
addpath ../sdwt2
addpath ../../CubeHelix/
addpath ../../export_fig_master

M = 1024;
N = 1024;
Qx = 4;
Qy = 4;
Q = 9; % number of spectral facets
tol = 1e-3;
d = 0.5*(M/Qx);

Q = Qx*Qy;
rg_y = domain_decomposition(Qy, M);
rg_x = domain_decomposition(Qx, N);
I = zeros(Q, 2);
dims = zeros(Q, 2);
for qx = 1:Qx
    for qy = 1:Qy
        q = (qx-1)*Qy+qy;
        I(q, :) = [rg_y(qy, 1)-1, rg_x(qx, 1)-1];
        dims(q, :) = [rg_y(qy,2)-rg_y(qy,1)+1, rg_x(qx,2)-rg_x(qx,1)+1];
    end
end

rg_yo = domain_decomposition_overlap2(Qy, M, d);
rg_xo = domain_decomposition_overlap2(Qx, N, d);
Io = zeros(Q, 2);
dims_o = zeros(Q, 2);
for qx = 1:Qx
    for qy = 1:Qy
        q = (qx-1)*Qy+qy;
        Io(q, :) = [rg_yo(qy, 1)-1, rg_xo(qx, 1)-1];
        dims_o(q, :) = [rg_yo(qy,2)-rg_yo(qy,1)+1, rg_xo(qx,2)-rg_xo(qx,1)+1];
    end
end

qx = 2;
qy = 2;
q = (qx-1)*Qy + qy;
tol = 1e-3;
if qx == 1
    wdx = [ones(1, dims(q,2)-d), linspace(1-tol, tol, d)];
elseif qx == Qx
    wdx = [linspace(tol, 1-tol, d), ones(1, dims_o(q,2)-2*d), ones(1, d)];
else
    wdx = [linspace(tol, 1-tol, d), ones(1, dims_o(q,2)-2*d), linspace(1-tol, tol, d)];
end    
if qy == 1
    wdy = [ones(1, dims(q,1)-d), linspace(1-tol, tol, d)];
elseif qy == Qy
    wdy = [linspace(tol, 1-tol, d), ones(1, dims_o(q,1)-2*d), ones(1, d)];
else
    wdy = [linspace(tol, 1-tol, d), ones(1, dims_o(q,1)-2*d), linspace(1-tol, tol, d)];
end
w = (wdy.').*wdx; 

%%

f = figure;
imagesc(w); colormap('cubehelix')
axis image;

% vertical lines
xx = [d+1, d+1];
yy = [1, dims_o(q,1)];
line(xx,yy,'Color','r','LineStyle','--', 'LineWidth', 2)

% horizontal lines
xx = [1, dims_o(q,2)];
yy = [d+1, d+1];
line(xx,yy,'Color','r','LineStyle','--', 'LineWidth', 2)

set(gca,'xticklabel',[]) % remove default ticks from imagesc
set(gca,'yticklabel',[])

mkdir figs
set(f,'units','normalized','outerposition',[0 0 0.52 0.9])
set(gca, 'Color', 'none'); % sets axes background
export_fig(['figs/weights.pdf'], '-transparent')
close(f)

f = figure;
imagesc(w); colormap('cubehelix')
axis image;
set(gca,'xticklabel',[]) % remove default ticks from imagesc
set(gca,'yticklabel',[])
set(f,'units','normalized','outerposition',[0 0 0.52 0.9])
set(gca, 'Color', 'none'); % sets axes background
export_fig(['figs/weights2.pdf'], '-transparent')
close(f)

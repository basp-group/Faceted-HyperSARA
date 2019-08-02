clc; clear all; close all;
format compact;

addpath ../../data
addpath ../../sdwt2
addpath ../../../subaxis
addpath ../../../export_fig_master

%% faceting parameters
N = [1024, 1024];
Qx = 4;
Qy = 4;
Q = Qx*Qy;
d = 100;

wlt_basis = {'db8'};
L = 16;
nlevel = 4;

% plot parameters
im_choice = 'W28';
bool_title = false;
clims_log = [-4, 0];

%% load image
x = fitsread(['W28_', num2str(N(1)),'.fits']);
x = flipud(x); % due to data reading in Matlab

%% define tessellation
% non-redundant image tessellation
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

% overlap for the sdwt2
[~, ~, I_overlap_ref, dims_overlap_ref, I_overlap, dims_overlap, ...
    ~, ~, status, offset, offsetL, offsetR, Ncoefs, temLIdxs, temRIdxs] = generate_segdwt_indices(N, I, dims, nlevel, wlt_basis, L);

% overlap for the nuclear norm
rg_yo = domain_decomposition_overlap2(Qy, N(1), d);
rg_xo = domain_decomposition_overlap2(Qx, N(2), d);
Io = zeros(Q, 2);
dims_o = zeros(Q, 2);
overlap = zeros(Q, 2);
for qx = 1:Qx
    for qy = 1:Qy
        q = (qx-1)*Qy+qy;
        Io(q, :) = [rg_yo(qy, 1)-1, rg_xo(qx, 1)-1];
        dims_o(q, :) = [rg_yo(qy,2)-rg_yo(qy,1)+1, rg_xo(qx,2)-rg_xo(qx,1)+1];
        overlap(q, :) =  dims_o(q,:) - dims(q,:);%max(dims_overlap_ref(q, :), dims_o(q, :)) - dims(q,:); % issue here! max(dims_overlap{q}, [], 1)
    end
end

%% plot the results
% figure;
% x_final = zeros(N);
% for q = 1:Q
%     [qy, qx] = ind2sub([Qy, Qx], q);
%     qm = (qy-1)*Qx + qx;
%     subplot(Qy, Qx, qm);
%     imagesc(log10(x(I(q,1)+1:I(q,1)+dims(q,1), I(q,2)+1:I(q,2)+dims(q,2))), [-5, 0]); 
%     %colormap gray; % to be checked again (problem in the way the images are displayed)
%     set(gca,'xticklabel',[]) % remove default ticks from imagesc
%     set(gca,'yticklabel',[])
%     axis image
% end
% set(gca, 'Color', 'none');
% export_fig(['figs/x_', name_img, '.pdf'], '-transparent')


%  nRows = 4;
%  nCols = 4;
%  % - Create figure, set position/size to almost full screen.
%  figure() ;
%  set( gcf, 'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.8] ) ;
%  % - Create grid of axes.
%  [blx, bly] = meshgrid( 0.05:0.9/nCols:0.9, 0.05:0.9/nRows:0.9 ) ;
%  hAxes = arrayfun( @(x,y) axes( 'Position', [x, y, 0.9*0.9/nCols, 0.9*0.9/nRows] ), blx, bly, 'UniformOutput', false ) ;
%  % - "Plot data", update axes parameters.
%  for q = 1 : numel( hAxes )
%     axes( hAxes{q} ) ;
%     [qy, qx] = ind2sub([Qy, Qx], q);
%     qm = (qy-1)*Qx + qx;
%     imagesc(log10(x(I(qm,1)+1:I(qm,1)+dims(qm,1), I(qm,2)+1:I(qm,2)+dims(qm,2))), [-5, 0]); 
%     %colormap gray; % to be checked again (problem in the way the images are displayed)
%     set(gca,'xticklabel',[]) % remove default ticks from imagesc
%     set(gca,'yticklabel',[])
%     axis image
%     set( gca, 'Visible', 'off' ) ;
%  end

%=========================================================================%
% Plot parameters
%=========================================================================%
fontsize=14;
orient = 'portrait';
nc = Qx; % nb of columns
nr = Qy; % nb of rows (number of methods)

% setting the Matlab figure
f=figure('visible','on');
set(gca, 'Color', 'none'); % Sets axes background
set(f,'PaperUnits','normalized')
set(f,'PaperType','A4');
set(f,'PaperOrientation',orient);
set(f,'units','normalized','outerposition',[0 0 0.40 0.70]) % 0.7 0.59
% colormap hot
% figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape');

for i=1:nr
    for ii=1:nc
        %ax = subaxis(nr,nc,ii,i,'SpacingHoriz',0.01,'SpacingVert',0.01); % best option !
        q = (ii-1)*nr + i;
        
        qm = (i-1)*nc + ii;
        subplot(nr, nc, qm)
%         qm = (qy-1)*Qx + qx;
        imagesc(log10(x(I(q,1)+1:I(q,1)+dims(q,1), I(q,2)+1:I(q,2)+dims(q,2))), [-5, 0]); 
        
%         if i == nr
%             set(l ,'Fontsize',fontsize)
%             pos = get(l, 'Position');
%             set(l, 'Position', pos + [0.0, 0.1, 0]);
%         end
        set(gca,'xticklabel',[]) % remove default ticks from imagesc
        set(gca,'yticklabel',[])
        axis image
        if bool_title && i==1
            if strcmp(im_choice, 'W28')
                title('','Interpreter','latex', 'Fontsize',fontsize)
            end
        end
        
%         if ii==nc
%             % colormap gray
%             % Get positions of all the subplot
%             posa = get(ax,'position');
%             h = colorbar;
%             %h.TicksMode = 'manual';
%             %h.TickLabels = {'10^-4','10^-2','1'};
%             %h.TickLabelInterpreter = 'latex';
%             set(h,'Fontsize',fontsize)
%             % Reset ax(3) position to before colorbar
%             set(ax,'position',posa)
%             % Set everything to units pixels (avoids dynamic reposition)
%             %                 set([ax h],'units','pix')
%             % Widen figure by a factor of 1.1 (tweak it for needs)
%             posf = get(gcf,'position');
%             set(gcf,'position',[posa(1:3) posf(4)])
%         end
    end
end

% setting the Matlab figure
f=figure('visible','on');
set(gca, 'Color', 'none'); % Sets axes background
set(f,'PaperUnits','normalized')
set(f,'PaperType','A4');
set(f,'PaperOrientation',orient);
set(f,'units','normalized','outerposition',[0 0 0.40 0.70])

total_size_x = sum(dims_o(:,2));
total_size_y = sum(dims_o(:,1));
max_dims_o = max(dims_o, [], 1);
min_dims_o = min(dims_o, [], 1);

% see how to properly display the images with their respective size (respect local and global ratio)

for i=1:nr
    for ii=1:nc
        if ii == 1
            pos_prev = [0, 0, 0, 0];
        end
        ax = subaxis(nr,nc,ii,i,'SpacingHoriz',0.01,'SpacingVert',0.01); % best option !
        % dims_o(q,2)/max_dims_o(2), dims_o(q,1)/max_dims_o(1)
        % min_dims_o(2)/dims_o(q,2), min_dims_o(1)/dims_o(q,1)
        q = (ii-1)*nr + i;
%         pos = get(ax, 'Position');
%         set(gca, 'Position', [pos(1)-pos_prev(1)+pos_prev(3), pos(2)-pos_prev(2)+pos_prev(4), dims_o(q,2)/(Qx*max_dims_o(2)), dims_o(q,1)/(Qy*max_dims_o(1))]) % ok for the dims, not ok for the positions
%         pos_prev = pos;
        
%         qm = (i-1)*nc + ii;
%         subplot(nr, nc, qm)
        imagesc(log10(x(Io(q,1)+1:Io(q,1)+dims_o(q,1), Io(q,2)+1:Io(q,2)+dims_o(q,2))), [-5, 0]); 
        hold on
        
        % vertical lines
        if overlap(q,2) == 0 && overlap(q,1) > 0
            xx = [overlap(q,2)+1, overlap(q,2)+1];
            yy = [1, dims_o(q,1)];
            line(xx,yy,'Color','r','LineStyle','--')
        else
            xx = [overlap(q,2), overlap(q,2)];
            yy = [1, dims_o(q,1)];
            line(xx,yy,'Color','r','LineStyle','--')
        end
        
        % horizontal lines
        if overlap(q,1) == 0 && overlap(q,2) > 0
            xx = [1, dims_o(q,2)];
            yy = [overlap(q,1)+1, overlap(q,1)+1];
            line(xx,yy,'Color','r','LineStyle','--')
        else
            xx = [1, dims_o(q,2)];
            yy = [overlap(q,1), overlap(q,1)];
            line(xx,yy,'Color','r','LineStyle','--')
        end

        set(gca,'xticklabel',[]) % remove default ticks from imagesc
        set(gca,'yticklabel',[])
        %axis image
        if bool_title && i==1
            if strcmp(im_choice, 'W28')
                title('','Interpreter','latex', 'Fontsize',fontsize)
            end
        end
    end
end

function [L,flag_ind,xinp_tot,xinp_multscale]  = create_inpainting_operator_test2_3D(Mask, Size_Gauss_kern,Size_Gauss_kern_freq, xmap_full, ch)

%% Initialization

Nx = size(Mask,2) ; Ny = size(Mask,1) ;

% Define the variance for the Gaussian weights
Sigma_Gauss_kern = 0.5* sqrt(Size_Gauss_kern) ;  %%%%%%%%%%%%%%%%%%%%%%%%%% This needs to be tuned

% save the inpainted images for each scale
xinp_multscale = cell(1,1,length(Size_Gauss_kern)) ;

%% create matrix

Nm = sum(Mask(:)) ; % dim variable inside mask
Nc = Nx*Ny - Nm ;   % dim variable outside mask

s = 0 ;

%%
for Sgauss = Size_Gauss_kern
s = s+1 ;
% build the Gaussian weights
sigmgauss = Sigma_Gauss_kern(s) ; 

% 2D
% wG = fspecial('gaussian',Sgauss*[1,1],sigmgauss);
% wG( floor(Sgauss/2)+1, floor(Sgauss/2)+1) = 0 ;
% wG = wG/( wG( floor(Sgauss/2), floor(Sgauss/2)+1)) ; %%%%%%% so that (+) neighbors are equal to one

% 3D
AA = zeros(Sgauss,Sgauss,Size_Gauss_kern_freq);
AA(round(Sgauss/2),round(Sgauss/2),round(Size_Gauss_kern_freq/2)) = 1;
wG = imgaussfilt3(AA,sigmgauss,'FilterSize',[Sgauss, Sgauss, Size_Gauss_kern_freq]);
wG(round(Sgauss/2),round(Sgauss/2),:) = 0 ;
ww = ones(1,1,3); % * 0.1;
ww(:,:,2) = 1;
wG = wG .* ww;
wG = wG/( wG( floor(Sgauss/2), floor(Sgauss/2)+1, floor(Size_Gauss_kern_freq/2)+1)) ; %%%%%%% so that (+) neighbors are equal to one

% special case for first and last channels
if ch - floor(Size_Gauss_kern_freq/2) < 1
wG(:,:,1:abs(ch - round(Size_Gauss_kern_freq/2))) = [] ;
end
if ch + floor(Size_Gauss_kern_freq/2) > size(xmap_full,3)
wG(:,:, end-abs(ch + floor(Size_Gauss_kern_freq/2)-size(xmap_full,3))+1:end) = [] ;
end

%% check intensity difference between channels (filter size 3 in frequency)
xmap_crop = xmap_full(:,:,max(ch - floor(Size_Gauss_kern_freq/2), 1) : min(ch + floor(Size_Gauss_kern_freq/2), size(xmap_full,3)));
ww_per = 0.1;

if size(xmap_crop,3) == Size_Gauss_kern_freq
    mask_3d = repmat(Mask,1,1,Size_Gauss_kern_freq);
    ss = mask_3d .* xmap_crop;
    
    for i = 1 : size(xmap_crop,3)
        tempo = ss(:,:,i);
        mm(i) = mean(tempo(:));
    end
    
    if abs(mm(2) - mm(1))/mm(2) >= ww_per && abs(mm(2) - mm(3))/mm(2) < ww_per
        wG(:,:,1) = [];
        xmap_crop(:,:,1) = [];
        flag_ind = [ch ch+1];
    elseif abs(mm(2) - mm(3))/mm(2) >= ww_per && abs(mm(2) - mm(1))/mm(2) < ww_per
        wG(:,:,3) = [];
        xmap_crop(:,:,3) = [];
        flag_ind = [ch-1 ch];
    elseif abs(mm(2) - mm(3))/mm(2) >= ww_per && abs(mm(2) - mm(1))/mm(2) >= ww_per
        wG = wG(:,:,2);
        xmap_crop = xmap_crop(:,:,2);
        flag_ind = ch;
    else
        flag_ind = [ch-1 ch ch+1];
    end
    
elseif ch == 1
    mask_3d = repmat(Mask,1,1,2);
    ss = mask_3d .* xmap_crop;
    
    for i = 1 : size(xmap_crop,3)
        tempo = ss(:,:,i);
        mm(i) = mean(tempo(:));
    end
    
    if abs(mm(2) - mm(1))/mm(2) >= ww_per 
        wG = wG(:,:,1);
        xmap_crop = xmap_crop(:,:,1);
        flag_ind = ch;
    else
        flag_ind = [ch ch+1];
    end
    
elseif ch == size(xmap_full,3)
    mask_3d = repmat(Mask,1,1,2);
    ss = mask_3d .* xmap_crop;
    
    for i = 1 : size(xmap_crop,3)
        tempo = ss(:,:,i);
        mm(i) = mean(tempo(:));
    end
    
    if abs(mm(2) - mm(1))/mm(2) >= ww_per
        wG = wG(:,:,2);
        xmap_crop = xmap_crop(:,:,2);
        flag_ind = ch;
    else
        flag_ind = [ch-1 ch];
    end
    
end

%wG = wG(:,:,3);
Sgauss_new = size(wG,3);


Maskt = Mask;
xinp = xmap_full(:,:,ch);
xinp(Mask>0)=0;
ro = 0 ;
indx = [] ; indy = [] ;
Sbord = 0 ;
% edges of Mask at last iteration
indx_old = indx ;
indy_old = indy ;
Sbord_old = Sbord ;

for dd = 1 : Sgauss_new 
Ltmp{dd} = sparse(Nm, Nc);
L{dd} = sparse(Nm, Nc);
end





%%
while sum(Maskt(:))>0
ro = ro+1 ;


%% mask update

% search indices in the Mask
indx_ = [] ;
indy_ = [] ;
for i = 1:Nx
for j = 1:Ny
if Maskt(j,i) ==1
indx_= [indx_; i] ;
indy_= [indy_; j] ;
end
end
end
Ind = [indy_,indx_];

% search the edges of the mask
Maskx = Maskt ; Maskx(2:end, 2:end) = Maskt(2:end,2:end) - Maskt(1:end-1,2:end) ;
for i=1:Nx
for j=1:Ny
if Maskx(j,i) == -1
Maskx(j,i) = 0 ;
Maskx(j-1,i) = 1 ;
end
end
end
Masky = Maskt ; Masky(2:end, 2:end) = Maskt(2:end,2:end) - Maskt(2:end,1:end-1) ;
for i=1:Nx
for j=1:Ny
if Masky(j,i) == -1
Masky(j,i) = 0 ;
Masky(j,i-1) = 1 ;
end
end
end
Maskg = min(abs(Maskx)+abs(Masky),1);
Sbord = sum(Maskg(:)) ; % number of elements in the edges

% search indices of the edges of the Mask
indx = [] ;
indy = [] ;
for i = 1:Nx
for j = 1:Ny
if Maskg(j,i) ==1
indx= [indx, i] ;
indy= [indy, j] ;
end
end
end


% 3D
xtmp = xmap_crop;
for i = 1 : size(xtmp,3)
    tmp = xtmp(:,:,i);
    tmp(Maskt>0) = 0;
    xtmp(:,:,i) = tmp;
end

% 2D
% xtmp = xmap_full(:,:,ch);
% xtmp(Maskt>0) = 0;

%%
for t=1:Sbord %each pixel in the edges of the mask
    
i=indx(t);
j=indy(t);

%%
% truncated neighborhood if the pixel of interest is too close to the edge 
% of the global image
Voisx = max(i - floor(Sgauss/2), 1) : min(i + floor(Sgauss/2), Nx) ;
Voisy = max(j - floor(Sgauss/2), 1) : min(j + floor(Sgauss/2), Ny) ;
Vois = [] ;
for u = 1:length(Voisx)
Vois = [Vois ; [Voisy', Voisx(u)*ones(length(Voisy),1)] ] ; %all Voisx, Voisy combinations
end

% keep neighboors only outside the Mask
ind_rem = ismember(Vois, Ind, 'rows') ;
Vois(ind_rem,:) = [] ;

%%
% truncated weights if the pixel of interest is too close to the edge of
% the global image
wGtmp = wG;
if i - floor(Sgauss/2)<1
wGtmp(:,1:abs(i - floor(Sgauss/2)+1),:) = [] ;
end
if j - floor(Sgauss/2)<1
wGtmp(1:abs(j - floor(Sgauss/2)+1),:,:) = [] ;
end
if i + floor(Sgauss/2)>Nx
wGtmp(:,end-abs(i + floor(Sgauss/2)-Nx)+1:end,:) = [] ;
end
if j + floor(Sgauss/2)>Ny
wGtmp(end-abs(j + floor(Sgauss/2)-Ny)+1:end,:,:) = [] ;
end

WGt = zeros([size(Mask) Sgauss_new]);
WGt(Voisy,Voisx,:) = wGtmp/sum(wGtmp(:));

WG = zeros([size(Mask) Sgauss_new]);
for u = 1:numel(Vois(:,1)) 
WG(Vois(u,1), Vois(u,2),:) = WGt(Vois(u,1), Vois(u,2),:);
end

WG = WG / sum(WG(:));

%% remove coefficients from WG within the mask
%  (write them w.r.t. the other coefficients)

[S_up, iv, isup] = intersect(Vois, [indy_old', indx_old'], 'rows', 'legacy') ;

for tt = 1:size(S_up,1)    % go through all the neighbors within the mask already treated
    ii = S_up(tt,2) ; %indx_old(isup(tt)) ;
    jj = S_up(tt,1) ; %indy_old(isup(tt)) ;
    
    % find the corresponding pixel n in L(n,:)
    tmp = zeros(Ny,Nx) ;
    tmp(jj,ii) = 1 ;
    tmp = tmp(Mask>0) ;
    tmp2 = 1:size(tmp,1) ;
    n = tmp2(tmp>0); %the index of the pixel
    
    for kk = 1 : size(WG,3)
    % 1. save the weight
    wt = WG(jj,ii,kk);
    
    % 2. get the pixel neighbors
    WGt = zeros(Ny*Nx,1);
    WGt(Mask==0) = Ltmp{kk}(n,:)/sum(Ltmp{kk}(n,:));
    WGt = reshape(WGt,Ny,Nx) ;
    
    % 3. update WG
    WG(jj,ii,kk) = 0; 
    WG(:,:,kk) = WG(:,:,kk) + wt * WGt; 
    end
    
    %%%%%%% summary: for every pixel in the mask already treated:
    % 1. we save its weight in the current kernel; w
    % 3. we add the outside neighbors that constructed this pixel weighted
    % with w 
end

values = WG.*xtmp;
xinp(j,i) = sum(values(:));


%% build matrix for coeff (i,j)

WG_v = reshape(WG(:),numel(WG(:))/Sgauss_new,Sgauss_new);
WGt = WG_v((Mask==0),:);

% index correspondance within the mask
tmp = zeros(Ny,Nx) ;
tmp(j,i) = 1 ;
tmp = tmp(Mask>0) ;
tmp2 = 1:size(tmp,1) ;
n = tmp2(tmp>0) ;

for dd = 1 : Sgauss_new 
Ltmp{dd}(n,:) = sparse(WGt(:,dd))' ;
end

end

%%
% update mask
Maskt = Maskt - Maskg ;
% edges of Mask at last iteration
indx_old = [indx_old, indx] ;
indy_old = [indy_old, indy] ;
Sbord_old = Sbord_old+Sbord ;

end

%% update multiscaled inpainted images 
xinp_multscale{s} = xinp;
for dd = 1 : Sgauss_new 
L{dd} = L{dd}+Ltmp{dd};
end

end
%% global inpainted image
xinp_tot_ = (1/length(Size_Gauss_kern)) * squeeze(sum(cell2mat(xinp_multscale),3)) ;
xinp_tot = xmap_full(:,:,ch);
xinp_tot(Mask>0) = xinp_tot_(Mask>0);

% for dd = 1 : Sgauss_new 
% L{dd} = L{dd}/length(Size_Gauss_kern) ;
% end

end
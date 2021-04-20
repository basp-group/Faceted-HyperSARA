function [x0,X0] = Generate_cube_W28(file_name,f,emission_lines)

disp '/*----------------------------------------------------------------*/'
disp ' example: Simulation of wide band radio interferometric data NEW'
disp '/*----------------------------------------------------------------*/'

%% Input model image
im=fitsread(file_name);
[Ny,Nx] = size(im);
% im=((imresize(im,[Ny,Nx],'nearest')));
im=(im+abs(im))./2;% Enforce positivity
im=im./max(im(:));
%imagesc(log10(flip(im)));

%% Generate sources
col = [22  3   192 180 191 218 148 173 81  141 104] * Nx/256;
row = [154 128 97  77  183 174 207 188 210 139 137] * Nx/256;
r =   [18  47  30  53  21  48  15  40  17  35  71]  * Nx/256;

n_src = length(col);
S = zeros(Ny*Nx,n_src+1);
imageSize = size(im);

foreground = zeros(size(im));
%count = 1;
%figure(1);
for i = 1 : n_src
    [xx,yy] = ndgrid((1:imageSize(1))-row(i),(1:imageSize(2))-col(i));
    mask = double((xx.^2 + yy.^2)<r(i)^2);
    %mask = (mask * count);
    %count = count +1;
    %imshow(mask); hold on;
    if i == 9
        mask9 = ~mask;
    end
    if i == 11
        mask10 = mask;
        mask = mask9 & mask10;
    end
    s = im.*mask;
    s(s<5e-3) = 0;
    S(:,i) = s(:);
    %fitswrite(s,['s' num2str(i) '.fits']);
    foreground = foreground + s;
    %imagesc(log10(s));
    %waitforbuttonpress;
end

%% Background source
s_background = im - foreground;
s_background(s_background<0) = 0;
S(:,end) = s_background(:);
%fitswrite(s_background,['s' num2str(i+1) '.fits']);

n_src = n_src + 1
%%

% figure(1),imagesc(log10(max(im,1e-3)));cb = colorbar; axis image; axis off; colormap jet; set(cb, 'fontsize', 22);
% set(gcf, 'Position', [0 400 800 800]);
% caxis([-3, 0]);

%% Build spectra

% alpha = [ -5 -2 0.5 0.3 -4 -3 0.1 0.5 -6 -7 0.5 -4.5];
% beta =  [ -5 -2  3   5.1  -4 -3  4   4  -6 -7  5  -4.5];

% alpha = [0.3 0.3 0.5 0.1 0.3 0.5 0.1 0.3 0.5 0.1 -10 0.5];
% beta =  [6 6 6 5 5 5 4 4 4 3 -10 5];

%alpha = [0.3 0.3 0.5 0.5 0.7 0.7 0.9 0.9 0.9 0.5 0.5 0.7];
%beta =  [3 3 3 3 3 3 3 3 3 2 2 2];

%alpha = [0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.5 0.5 0.7];
%beta =  [3 3 3 3 3 3 3 3 3 2 2 2];

alpha = 0.2;
beta = 2;

alpha = repmat(alpha,1,n_src);
beta = repmat(beta,1,n_src);

%! test P.-A.
% alpha = [ -3 -2.5 0.5 0.3 -2 -1 0.1 0.5 -0.25 -4 0.5 2.5];
% beta =  [ -3 -2  3   5.1  -3 -2.5  4   4  0.25 -3  3  3];

c = length(f);
%
H = zeros(c,n_src);
HI = ones(c,n_src);

if emission_lines
    k = 3;
    for i=1:n_src
        HI(k,i) = 2;
        HI(c-k,i) = 1.5;
        HI(c-k-1,i) = 2;
        HI(c-k-2,i) = 1.5;
        HI(c-k-3,i) = 2;
        HI(c-k-4,i) = 1.5;
        k = k+3;
    end
end

police={'-*b','-*b','-*b'}; %, '-or', '-dg', '-sc','.-m', '-+c', '-hb', '-dr','*g', 'sk', ' -.-y' , '-*b'};

%%

%figure(2), hold on;

for i=1:length(alpha)
    
    H(:,i) = (f./f(1)).^(- alpha(i) + (beta(i).*log(f./f(1))));
    if emission_lines
        H(:,i) = H(:,i) .* HI(:,i);
    end
    %subplot 121; imagesc(log10(reshape(S(:,i),[1024 1024])));
    %subplot 122; plot(H(:,i),police{2});hold on; xlim([1 c]); ylim([1 6]);
    %     if i == 11
    %         plot(H(:,i),police{2}); ylim([1 25]);
    %     elseif i == 12
    %         plot(H(:,i),police{3}); ylim([1 10]);
    %     elseif i == 1 || i == 2 || i == 3
    %         plot(H(:,i),police{1}); ylim([1 25]);
    %     else
    %         plot(H(:,i),police{1}); ylim([1 10]);
    %     end
    %waitforbuttonpress;
    
    
end
%hold off;

%H(:,1) = H(:,1) - 1;

%% Simulation of data cube

X0 = S*H';
X0(X0<0) = 0;
% X0 = X0./max(X0(:));

for i=1:size(X0,2)
    x0(:,:,i)=reshape(X0(:,i),[Ny,Nx]);
end

rank(X0)

% figure(3),imagesc(log10(max(X0,1e-3)));cb = colorbar; colormap jet; set(cb, 'fontsize', 22);
% set(gcf, 'Position', [0 400 800 800]);
% caxis([-3, 0]);

%%
% figure(4),hold on
% for i = 1 : 10 : size(X0,1)
%     if X0(i,end) > 0
%         plot (X0(i,:))
%     end
% end
% xlim([1 c]);

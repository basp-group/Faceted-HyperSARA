clear;  close all
%% Freq info
nu_1 = 2.052e9; %starting freq
dnu = 16e6;     % freq step
L = 100   ;     %number of channels
nu_vect =[nu_1 (dnu*(1:L-1)+nu_1)];
%% for info
dfreq =5
nu_2k = nu_vect(1:dfreq:end);
nu_1k = nu_vect(1:dfreq:end);
nu_512 = nu_vect;


%% Cyg A images at different resolutions
im_ref = fitsread('cygA_Sband_Dabbech2021.fits');
Nh =3000; 
imr = imresize(im_ref,[Nh,Nh]);
im_2k= imr(Nh/2 +1-512:Nh/2 +1 +511,Nh/2 +50 -1024:Nh/2 +50 +1023);
im_2k = im_2k .*(im_2k>2e-5);
fitswrite(im_2k,['cygASband_' num2str(size(im_2k,1)) '_' num2str(size(im_2k,2)) '.fits'])
im_1k =imresize(im_2k,[512,1024]);
im_1k=im_1k.*(im_1k>4e-5);
fitswrite(im_1k,['cygASband_' num2str(size(im_1k,1)) '_' num2str(size(im_1k,2)) '.fits'])

im_512 =imresize(im_2k,[256,512]);
im_512=im_512.*(im_512>8e-5);
fitswrite(im_512,['cygASband_' num2str(size(im_1k,1)) '_' num2str(size(im_1k,2)) '.fits'])

% % figure,subplot 311, imagesc(log10(im_2k)),axis image ,colorbar, title ('Cyg A : 2kx1k')
% % subplot 312, imagesc(log10(im_1k)),axis image, colorbar, title ('Cyg A : 1kx512')
% % subplot 313, imagesc(log10(im_512)),axis image,colorbar, title ('Cyg A : 512x256')

im_512_init = im_512;
im_1k_init = im_1k;
im_2k_init = im_2k;
%
fact_dim = size(im_2k,1)./[size(im_2k,1) size(im_1k,1) size(im_512,1)];
%% spectral index map: 1st component:smooth image and then take the log

N_ =[64,64];
[X,Y] = meshgrid(-N_(1)/2:N_(1)/2 -1, -N_(1)/2:N_(2)/2  -1);
h= [X(:) Y(:)];
theta=0;
% 2kx1k
sig_x0=10;
sig2_y0=10;
a =cos(theta).^2/(2*sig_x0.^2) +sin(theta).^2/(2*sig2_y0.^2);
b=sin(2*theta).^2/(4*sig_x0.^2) +sin(2*theta).^2/(4*sig2_y0.^2);
c =sin(theta).^2/(2*sig_x0.^2) +cos(theta).^2/(2*sig2_y0.^2);
rho = @(h) exp(-(a*h(:,1).^2 +2*b.*h(:,1).*h(:,2)+c*h(:,2).^2));
beam_2k = reshape(rho(h),N_(1),N_(2));
beam_2k = beam_2k./sum(sum(beam_2k));
fwhm1 = 2*sqrt(2*log(2))*sig_x0;
fwhm2 = 2*sqrt(2*log(2))*sig2_y0;
im_2kSLog = log10(conv2(im_2k,beam_2k,'same').*(im_2k_init>0));
%1kx512
sig_x=sig_x0/fact_dim(2);
sig_y=sig2_y0/fact_dim(2);
a = cos(theta).^2/(2*sig_x.^2) +sin(theta).^2/(2*sig_y.^2);
b = sin(2*theta).^2/(4*sig_x.^2) +sin(2*theta).^2/(4*sig_y.^2);
c = sin(theta).^2/(2*sig_x.^2) +cos(theta).^2/(2*sig_y.^2);
rho = @(h) exp(-(a*h(:,1).^2 +2*b.*h(:,1).*h(:,2)+c*h(:,2).^2));
beam_1k = reshape(rho(h),N_(1),N_(2));
beam_1k = beam_1k./sum(sum(beam_1k));
im_1kSLog = log10(conv2(im_1k,beam_1k,'same').*(im_1k_init>0));
%512x256
sig_x = sig_x0/fact_dim(3);
sig_y = sig2_y0/fact_dim(3);
a = cos(theta).^2/(2*sig_x.^2) +sin(theta).^2/(2*sig_y.^2);
b = sin(2*theta).^2/(4*sig_x.^2) +sin(2*theta).^2/(4*sig_y.^2);
c = sin(theta).^2/(2*sig_x.^2) +cos(theta).^2/(2*sig_y.^2);
rho = @(h) exp(-(a*h(:,1).^2 +2*b.*h(:,1).*h(:,2)+c*h(:,2).^2));
beam_512 = reshape(rho(h),N_(1),N_(2));
beam_512 = beam_512./sum(sum(beam_512));
im_512SLog = log10(conv2(im_512,beam_512,'same').*(im_512_init>0));


%% spectral index map: 2nd component : gaussian process

sig_gprocess =10;%% sigma value for Gaussian process
N2k = size(im_2k);
theta=max(min(randn(1,1),1),-1)*pi;
sig_x=sig_gprocess/fact_dim(1);
sig_y=sig_gprocess/fact_dim(1);
a =cos(theta).^2/(2*sig_x^2) +sin(theta).^2/(2*sig_y^2);
b=sin(2*theta).^2/(4*sig_x^2) +sin(2*theta).^2/(4*sig_y^2);
c =sin(theta).^2/(2*sig_x^2) +cos(theta).^2/(2*sig_y^2);
rho = @(h) exp(-(a*h(1)^2 +2*b*h(1)*h(2)+c*h(2)^2));
[f1Index_2k,~]=stationary_Gaussian_process(N2k(1),N2k(2),rho);
f1Index_2k =f1Index_2k./max(max(abs(f1Index_2k)))  .*(im_2k_init>0);
%
N1k = size(im_1k);
sig_x=sig_gprocess/fact_dim(2);
sig_y=sig_gprocess/fact_dim(2);
a =cos(theta).^2/(2*sig_x^2) +sin(theta).^2/(2*sig_y^2);
b =sin(2*theta).^2/(4*sig_x^2) +sin(2*theta).^2/(4*sig_y^2);
c =sin(theta).^2/(2*sig_x^2) +cos(theta).^2/(2*sig_y^2);
rho = @(h) exp(-(a*h(1)^2 +2*b*h(1)*h(2)+c*h(2)^2));
[f1Index_1k,~]=stationary_Gaussian_process(N1k(1),N1k(2),rho);
f1Index_1k =f1Index_1k./max(max(abs(f1Index_1k)))  .*(im_1k_init>0);
%
N512= size(im_512);
sig_x=sig_gprocess/fact_dim(3);
sig_y=sig_gprocess/fact_dim(3);
a =cos(theta).^2/(2*sig_x^2) +sin(theta).^2/(2*sig_y^2);
b=sin(2*theta).^2/(4*sig_x^2) +sin(2*theta).^2/(4*sig_y^2);
c =sin(theta).^2/(2*sig_x^2) +cos(theta).^2/(2*sig_y^2);
rho = @(h) exp(-(a*h(1)^2 +2*b*h(1)*h(2)+c*h(2)^2));
[f1Index_512,~]=stationary_Gaussian_process(N512(1),N512(2),rho);
f1Index_512 =f1Index_512./max(max(abs(f1Index_512)))  .*(im_512_init>0);

% spectral index map
scale = 1;
spectralIds2k  = (scale .*f1Index_2k.* (im_2k_init>0)+im_2kSLog);
spectralIds1k  = (scale .*f1Index_1k.* (im_1k_init>0)+im_1kSLog);
spectralIds512  = (scale .*f1Index_512.* (im_512_init>0)+im_512SLog);

spectralIds2k =(spectralIds2k -max(max(spectralIds2k)) -0.3);
spectralIds1k =(spectralIds1k -max(max(spectralIds1k)) -0.3);
spectralIds512 =(spectralIds512 -max(max(spectralIds512)) -0.3);
spectralIds2k(spectralIds2k<-5) =NaN;
spectralIds1k(spectralIds1k<-5) =NaN;
spectralIds512(spectralIds512<-5) =NaN;

figure,subplot 311, imagesc(spectralIds2k),axis image,colorbar,title('Spectral index')
subplot 312, imagesc(spectralIds1k),axis image,colorbar,title('Spectral index')
subplot 313, imagesc(spectralIds512),axis image,colorbar,title('Spectral index')


%% curvature map
scaleC =0.5;
theta= max(min(randn(1,1),1),-1)*pi;
sig_x0 = sig_gprocess /fact_dim(1);
sig_y = sig_gprocess /fact_dim(1);
a =cos(theta).^2/(2*sig_x0^2) +sin(theta).^2/(2*sig_y^2);
b = sin(2*theta).^2/(4*sig_x0^2) +sin(2*theta).^2/(4*sig_y^2);
c =sin(theta).^2/(2*sig_x0^2) +cos(theta).^2/(2*sig_y^2);
rho = @(h) exp(-(a*h(1)^2 +2*b*h(1)*h(2)+c*h(2)^2));
[f1Index_2k,~]=stationary_Gaussian_process(N2k(1),N2k(2),rho);
Curv2k = scaleC *f1Index_2k./max(max(abs(f1Index_2k))) .*(im_2k_init>0);
Curv2k(Curv2k<-scaleC) =NaN;
%
% theta= max(min(randn(1,1),1),-1)*pi;
sig_x0 = sig_gprocess /fact_dim(2);
sig_y = sig_gprocess /fact_dim(2);
a =cos(theta).^2/(2*sig_x0^2) +sin(theta).^2/(2*sig_y^2);
b = sin(2*theta).^2/(4*sig_x0^2) +sin(2*theta).^2/(4*sig_y^2);
c =sin(theta).^2/(2*sig_x0^2) +cos(theta).^2/(2*sig_y^2);
rho = @(h) exp(-(a*h(1)^2 +2*b*h(1)*h(2)+c*h(2)^2));
[f1Index_1k,~]=stationary_Gaussian_process(N1k(1),N1k(2),rho);
Curv1k = scaleC *f1Index_1k./max(max(abs(f1Index_1k))) .*(im_1k_init>0);
Curv1k(Curv1k<-scaleC) =NaN;
%
% theta= max(min(randn(1,1),1),-1)*pi;
sig_x0 = sig_gprocess /fact_dim(3);
sig_y = sig_gprocess /fact_dim(3);
a =cos(theta).^2/(2*sig_x0^2) +sin(theta).^2/(2*sig_y^2);
b = sin(2*theta).^2/(4*sig_x0^2) +sin(2*theta).^2/(4*sig_y^2);
c =sin(theta).^2/(2*sig_x0^2) +cos(theta).^2/(2*sig_y^2);
rho = @(h) exp(-(a*h(1)^2 +2*b*h(1)*h(2)+c*h(2)^2));
[f1Index_512,~]=stationary_Gaussian_process(N512(1),N512(2),rho);
Curv512 = scaleC *f1Index_512./max(max(abs(f1Index_512))) .*(im_512_init>0);
Curv512(Curv512<-scaleC) =NaN;

figure,subplot 311, imagesc(Curv2k),axis image,colorbar,title('Curv')
subplot 312, imagesc(Curv1k),axis image,colorbar,title('Curv')
subplot 313, imagesc(Curv512),axis image,colorbar,title('Curv')

%% Build cube

cube2k= zeros(N2k(1),N2k(2),L);
cube1k= zeros(N1k(1),N1k(2),L);
cube512= zeros(N512(1),N512(2),L);

spectralIds2k(isnan(spectralIds2k)) =0;
spectralIds1k(isnan(spectralIds1k)) =0;
spectralIds512(isnan(spectralIds512)) =0;

Curv2k(isnan(Curv2k)) =0;
Curv1k(isnan(Curv1k)) =0;
Curv512(isnan(Curv512)) =0;

for i =1:numel(nu_vect)
    cube2k(:,:,i) = im_2k_init.*(nu_vect(i)./nu_1).^(spectralIds2k+Curv2k*log(nu_vect(i)./nu_1) );
    cube1k(:,:,i) = im_1k_init.*(nu_vect(i)./nu_1).^(spectralIds1k+Curv1k*log(nu_vect(i)./nu_1) );
    cube512(:,:,i) = im_512_init.*(nu_vect(i)./nu_1).^(spectralIds512+Curv512*log(nu_vect(i)./nu_1) );

end
cube1k =cube1k(:,:,1:dfreq:end);
cube2k =cube2k(:,:,1:dfreq:end);

fitswrite(cube2k,['cygASband_Cube_' num2str(size(cube2k,1)) '_' num2str(size(cube2k,2)) '_' num2str(size(cube2k,3))   '.fits']);
fitswrite(cube1k,['cygASband_Cube_' num2str(size(cube1k,1)) '_' num2str(size(cube1k,2)) '_' num2str(size(cube1k,3))   '.fits']);
fitswrite(cube512,['cygASband_Cube_' num2str(size(cube512,1)) '_' num2str(size(cube512,2))  '_' num2str(size(cube512,3))   '.fits']);

fitswrite(spectralIds2k,['CurvMap' num2str(size(cube2k,1)) '_' num2str(size(cube2k,2)) '_' num2str(size(cube2k,3))     '.fits']);
fitswrite(spectralIds1k,['CurvMap' num2str(size(cube1k,1)) '_' num2str(size(cube1k,2))  '_' num2str(size(cube1k,3))    '.fits']);
fitswrite(spectralIds512,['CurvMap' num2str(size(cube512,1)) '_' num2str(size(cube512,2)) '_' num2str(size(cube512,3))     '.fits']);

fitswrite(Curv2k,['SpectralIdxMap' num2str(size(cube2k,1)) '_' num2str(size(cube2k,2))  '_' num2str(size(cube2k,3))    '.fits']);
fitswrite(Curv1k,['SpectralIdxMap' num2str(size(cube1k,1)) '_' num2str(size(cube1k,2)) '_' num2str(size(cube1k,3))     '.fits']);
fitswrite(Curv512,['SpectralIdxMap' num2str(size(cube512,1)) '_' num2str(size(cube512,2))  '_' num2str(size(cube512,3))    '.fits']);



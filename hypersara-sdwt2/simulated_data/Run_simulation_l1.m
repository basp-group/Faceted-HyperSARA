function [param_data, param_algo_map, result_MAP, param_hpd, param_inpaint, result] ...
          = Run_simulation_l1...
          (save_name, im_name, size_im, cov_type, Tinst, sig2,seed, nlevel, wlt_basis,...
          Cropx, Cropy, Mincrop, Maxcrop, Test_choice)

      
%% Load original image

if isempty(size_im)
[im, N, Ny, Nx] = util_read_image_(im_name);
elseif numel(size_im) == 1
[im, N, Ny, Nx] = util_read_image_(im_name, size_im, size_im);
elseif numel(size_im) == 2
[im, N, Ny, Nx] = util_read_image_(im_name, size_im(1), size_im(2));
else
disp('problem with size')
return;
end

param_data.im = im ;
param_data.N = N; % number of pixels in the image
param_data.Nx = Nx ;
param_data.Ny = Ny ;


%% Generate uv coverage -- continuous Fourier under-sampling --

rng(seed);

param_data.sigma2 = sig2 ;
param_data.sigma_noise = sqrt(sig2) ;

% % % % Generate Gaussian random u-v sampling
% % % param_data.p = 0.5 ;
% % % param_data.sigma = pi/4 ;
% % % [u, v] = util_gen_sampling_pattern('gaussian', param_data);

param_data.cov_type = cov_type ;
param_data.T = Tinst ;
param_data.na = 100 ; %in case of random distribution of antennas

% ----------------------------------------------------- %
switch param_data.cov_type                              %
% ----------------------------------------------------- %
% VLA coverage ---------------------------------------- %
% ----------------------------------------------------- %
case 'vlaa'                                             %
fileID = fopen('vlaa.itrf-copie.txt') ;                 %
save_cov_file = ['gen_uv_tracks/ant_vla_pos.mat'];      %
% ----------------------------------------------------- %
% ASKA coverage --------------------------------------- %
% ----------------------------------------------------- %
case 'askap'                                            %
fileID = fopen('askap.itrf-copie.txt') ;                %
save_cov_file = ['gen_uv_tracks/ant_askap_pos.mat'];    %
% ----------------------------------------------------- %
% MeerKAT coverage ------------------------------------ %
% ----------------------------------------------------- %
case 'meerkat'                                          %
fileID = fopen('MeerKAT.enu.txt') ;                     %
save_cov_file = ['gen_uv_tracks/ant_meerkat_pos.mat'];  %
% ----------------------------------------------------- %
% random cont. antenna positions                        %
% ----------------------------------------------------- %
case 'random'                                           %
uv(1,:) = rand(1,2) ;                                   %
for alpha = 2:param_data.na                             %
uv_ = rand(1,2) ;                                       %
while ismember(uv_,uv,'rows')                           %
uv_ = rand(1,2) ;                                       %
end                                                     %
uv(alpha,:) = uv_ ;                                     %
end                                                     %
Pos_ant = 1e06*[uv, zeros(param_data.na,1)] ;           %
save_cov_file = ['gen_uv_tracks/rand_pos.mat'];         %
% ----------------------------------------------------- %
end                                                     %
% ----------------------------------------------------- %
if strcmp(param_data.cov_type, 'rand_ant') == 0
C = textscan(fileID,'%f %f %f %s %s %s');
Pos_ant = cell2mat({C{1} C{2} C{3}}) ;
param_data.na = max(size(Pos_ant)) ;
fclose(fileID);
end

save(save_cov_file,'Pos_ant')

%%
param_data.M = param_data.T*param_data.na*(param_data.na-1)/2 ;
param_data.ant_pair = param_data.na*(param_data.na-1)/2 ;

disp('--------------------------------------------')
disp(['nombre de mesures : ',num2str(param_data.M)])
disp(['nombre d antennes : ',num2str(param_data.na)])
disp(['nombre de paires d antennes : ',num2str(param_data.ant_pair)])
disp(['nombre de pts / paires d antennes : ',num2str(param_data.T)])
disp('--------------------------------------------')


% -------------------------------------------------------------------------
h = linspace(-6, 6, param_data.T)*pi/12;% hour angle range of +/-  hours
% % % h = linspace(-3, 3, param_data.T)*pi/12;% hour angle range of +/-  hours
dec = (pi/180)*(40);  % Cas A is 56.4
lat = (pi/180)*(38. + 59./60. + 48./3600.);  % College Park
% position de reference x0 dans le plan x-y
x0 = [mean(Pos_ant(:,1)), mean(Pos_ant(:,2))] ;
% -------------------------------------------------------------------------
[U,V,~] = generate_uv_cov_antennas(Pos_ant,x0,h,lat,dec, param_data.T) ;

uab = cell(1,param_data.T) ;
vab = cell(1,param_data.T) ;
for t=1:param_data.T
    uab{t} = zeros(param_data.ant_pair, 1) ;
    vab{t} = zeros(param_data.ant_pair, 1) ;
    i = 1 ;
    for alp = 1:param_data.na
        for bet = alp+1:param_data.na
            uab{t}(i) = U{alp}(t)- U{bet}(t) ;
            vab{t}(i) = V{alp}(t)- V{bet}(t) ;
            i=i+1 ;
        end
    end
end

uu = cell2mat(uab) ; vv = cell2mat(vab) ;
maxuv = max(max(abs(uu(:))), max(abs(vv(:)))) ;
uuu = uu/maxuv*pi*0.95 ; vvv = vv/maxuv*pi*0.95 ;

%%
% u = cell2mat(U) ; u = u(:) ;
% v = cell2mat(V) ; v = v(:) ;
% maxuv = max(max(abs(u)), max(abs(v))) ;
% u = u/maxuv*pi ;
% v = v/maxuv*pi ;

u = uuu(:);
v = vvv(:) ;


% Initialize nuFFT operator
% Generate measurment operator with nufft
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 8; % number of neighbours for nufft
Ky = 8; % number of neighbours for nufft
[A, AT, Gw, ~] = op_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);

param_data.Phi =@(x) Gw * A(x) ;
param_data.Phit =@(y) AT(Gw' * y) ;

% norm of the measurement operator
param_data.normPhi = op_norm(param_data.Phi, param_data.Phit, [Ny, Nx], 1e-4, 200, 0);    


% Generate noisy measurements
y0 = param_data.Phi(im);
noise = param_data.sigma_noise*(randn(size(y0)) + 1i*randn(size(y0))) ;
param_data.y = y0 + noise;


param_data.M = length(y0) ;

% param_data.l2bound = 1.1 * norm(noise) ;
param_data.l2bound = sqrt(2*param_data.M + 2* sqrt(4*param_data.M)) *  param_data.sigma_noise ;

% Sparsity basis for l1 regularization
[Psi, Psit] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
normPsi = 1 ; %if normalized



%% Solve the MAP problem



param_algo_map.NbIt = 60000 ;
param_algo_map.stop_crit = 1e-6 ;
param_algo_map.Psit = Psit ;
param_algo_map.Psi = Psi ;
param_algo_map.normPsi = normPsi ;
param_algo_map.x0 = 0 * real(param_data.Phit(param_data.y)) ;
param_algo_map.x0 = max(param_algo_map.x0 / max(param_algo_map.x0(:)), 0) ;
param_algo_map.lambda = 1e-3 ; 
param_algo_map.gammat = param_data.normPhi/1e4 ; 
param_algo_map.display = 500 ;

tic
result_MAP = solve_MAP_constrainedPB(param_data,param_algo_map) ;
result_MAP.comp_time = toc ;

figure, 
subplot 221
imagesc(log10(result_MAP.x)), axis image; colorbar, colormap jet, caxis([-3.5,0])
xlabel('MAP log')
subplot 222
imagesc(log10(im)), axis image; colorbar, colormap jet, caxis([-3.5,0])
xlabel('true log')
subplot 223
imagesc(result_MAP.x), axis image; colorbar, colormap jet, caxis([0,1])
xlabel('MAP lin')
subplot 224
imagesc(im), axis image; colorbar, colormap jet, caxis([0,1])
xlabel('true lin')


save([save_name,'_MAP.mat'], ...
    'im_name', 'sig2','seed', 'nlevel', 'wlt_basis', ...
    'param_data', 'param_algo_map', 'result_MAP') ;

%% HPD definition and parameters

xmap = result_MAP.x ;

param_hpd.lambda_t = param_data.N / sum(abs(param_algo_map.Psit(xmap))) ;


alpha = 1e-2 ; 
talpha = sqrt( (16*log(3/alpha)) / param_data.N );
HPDconstraint = param_hpd.lambda_t* sum(abs(param_algo_map.Psit(xmap))) ...
                + param_data.N*(1+talpha);
param_data.HPDconstraint = HPDconstraint/param_hpd.lambda_t ;




param_hpd.NbIt = 3000 ;
param_hpd.Psit = param_algo_map.Psit ;
param_hpd.Psi = param_algo_map.Psi ;
param_hpd.normPsi = param_algo_map.normPsi ;
param_hpd.lambda = param_algo_map.lambda ;
param_hpd.display = 200 ;
param_hpd.perc = 1 ;
param_hpd.cond_stop = 1e-4 ;


%% Structure smoothing definition and parameters


result = cell(length(Test_choice),1) ;
for count_test = 1:length(Test_choice)
test_choice = Test_choice(count_test) ;
param_algo.test_choice = test_choice ;


cropx = Cropx(test_choice,1) : Cropx(test_choice,2) ;
cropy = Cropy(test_choice,1) : Cropy(test_choice,2) ;
mincrop = Mincrop(test_choice) ;
maxcrop = Maxcrop(test_choice) ;



if maxcrop==Inf
Mask = zeros(size(xmap));
Mask(cropy,cropx) = 1;
Mask = Mask.*(xmap>mincrop);
Mask = Mask.*(xmap<maxcrop);
Mask = imdilate(Mask,strel('disk',3));
else
Mask = ones(size(xmap)) ;
Mask(xmap<maxcrop) = 0 ;
Mask = imdilate(Mask,strel('disk',7));
Mask = abs(Mask-1);
end

tmp = xmap ;
tmp(Mask>0) = 0 ;
figure,
subplot 221, imagesc(xmap), axis image; colorbar, caxis([0,1]), colormap jet
subplot 222, imagesc(log10(xmap)), axis image; colorbar, caxis([-4.5 0]), colormap jet
subplot 223, imagesc(log10(tmp)), axis image; colorbar, caxis([-4.5 0]), colormap jet
subplot 224, imagesc(Mask), axis image; colorbar

  
param_inpaint.Mask = Mask ;

%%
if maxcrop<Inf
param_inpaint.choice_set = 'back0' ;

xmap_inp = xmap ;
xmap_inp(param_inpaint.Mask>0) = 0 ;
figure, imagesc(log10(xmap_inp)), axis image ; colorbar, colormap jet, caxis([-4.5 0])
else
param_inpaint.choice_set = 'l2_const' ;
param_inpaint.Size_Gauss_kern = [3,7,11] ;
disp('create inpainting operator...')
param_inpaint.L = create_inpainting_operator_test(param_inpaint.Mask, param_inpaint.Size_Gauss_kern, xmap) ;
disp('...done')
param_inpaint.L = sparse(param_inpaint.L) ;


param_inpaint.Lbar = sparse([-speye(sum(param_inpaint.Mask(:))), param_inpaint.L]) ;
param_inpaint.Nout = sum(param_inpaint.Mask(:)==0) ;
param_inpaint.normL = op_norm(@(x) param_inpaint.L*x, @(x) param_inpaint.L'*x, [param_inpaint.Nout,1], 1e-4, 200, 0);  
param_inpaint.normLbar = op_norm(@(x) param_inpaint.Lbar*x, @(x) param_inpaint.Lbar'*x, [numel(Mask),1], 1e-4, 200, 0);  

param_inpaint.Li =@(u) [param_inpaint.L*u ; u] ;
param_inpaint.Lit =@(x) param_inpaint.L'*x(param_inpaint.Mask>0) + x(param_inpaint.Mask==0) ;
param_inpaint.normLi = op_norm(param_inpaint.Lit, param_inpaint.Li, [numel(param_inpaint.Mask),1], 1e-4, 200, 0);  

xmap_inp = xmap ;
xmap_inp(param_inpaint.Mask>0) = param_inpaint.L*xmap(param_inpaint.Mask==0) ;
figure, imagesc(log10(xmap_inp)), axis image ; colorbar, colormap jet, caxis([-4.5 0])

param_inpaint.l2_mean = 0 ;
param_inpaint.l2_bound = sqrt(sum(abs(xmap_inp(param_inpaint.Mask>0)).^2)) ; % 1e-3 * sqrt(sum(param_inpaint.Mask(:))) ; 

Mop = sparse(sum(Mask(:)), numel(Mask)) ;
Mopcomp = sparse(numel(Mask)-sum(Mask(:)), numel(Mask)) ;
i = 1; ic = 1;
for n = 1:numel(Mask)
    if(Mask(n))==0
        Mopcomp(ic,n) = 1 ;
        ic = ic+1 ;
    else
        Mop(i,n) = 1;
        i = i+1 ;
    end
end

param_inpaint.Mop = Mop ;
param_inpaint.Mbarop = Mop - param_inpaint.L*Mopcomp;

xM_ = param_inpaint.Mop * xmap(:) ;
param_inpaint.l2b_Mask = sqrt(numel(xM_) + 2* sqrt(numel(xM_))) *  param_data.sigma_noise ;
end





param_inpaint.Psi = param_hpd.Psi ;
param_inpaint.Psit = param_hpd.Psit ;
param_inpaint.normPsi = param_hpd.normPsi ;
param_inpaint.l1bound = param_data.HPDconstraint ; 


param_inpaint.NbIt = 3000 ;
param_inpaint.display = 10 ;


%% POCS algorithm

param_algo.NbIt = 60000 ;
param_algo.stop_dist = 1e-5 ;
param_algo.stop_norm2 = 1e-5 ;
param_algo.stop_norm1 = 1e-5 ;
param_algo.stop_err_smooth = 1e-5 ;




l1_inpaint = sum(abs(param_algo_map.Psit(xmap_inp))) ;
l2_inpaint = sqrt(sum( abs( param_data.Phi(xmap_inp) - param_data.y ).^2 )) ;


disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp('*******************************************************')
disp('*******************************************************')
disp(['T = ',num2str(Tinst)])
disp(['sigma2 = ',num2str(param_data.sigma2)])
disp(['test choice no = ',num2str(test_choice)])
disp(['l1 inpaint = ',num2str(l1_inpaint)])
disp(['HPD bound = ',num2str(param_data.HPDconstraint)])
disp(['l2 data inpaint = ',num2str(l2_inpaint)])
disp(['data bound = ',num2str(param_data.l2bound)])


%%
if l1_inpaint <= param_data.HPDconstraint && l2_inpaint <= param_data.l2bound
    disp('inpainted image inside the HPD set')
    disp('*******************************************************')
    result{count_test}.xinp = xmap_inp ;
    result{count_test}.xhpd = xmap_inp ;
else
    disp('inpainted image OUTSIDE the HPD set')
    disp('*******************************************************')
    disp(' ')
% % %     result{count_test} = POCS_PD_global(xmap, xmap_inp, param_algo, param_data, param_hpd, param_inpaint) ;
    result{count_test} = POCS_PD_global_relax(xmap, xmap_inp, param_algo, param_data, param_hpd, param_inpaint) ;
end


   
rho  = norm( result{count_test}.xC(:) - result{count_test}.xS(:) ) / ...
       norm( xmap(:) - xmap_inp(:) ) ;
disp(' ')
disp('*****************************************')
disp(['T = ',num2str(Tinst), ', sig2 = ',num2str(sig2), ', seed = ',num2str(seed)])
disp(['test choice no. ',num2str(test_choice)])
disp(['       rho = ',num2str(rho)])
disp('*****************************************')


save([save_name,'_test=',num2str(test_choice),'.mat'], ...
    'im_name', 'sig2','seed', 'nlevel', 'wlt_basis','test_choice', 'Test_choice', 'count_test', ...
    'param_algo','param_data','param_hpd','param_inpaint','param_algo_map', ...
    'result_MAP', 'result', 'xmap_inp', 'rho') ;
    

end

end

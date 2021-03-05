function util_load_real_vis_data(nb_aggch)


speed = 299792458;
load(['CGCG044_046_FULL.mat']);

%% Split data
% number of overall channels
c = size(vis,2);
fprintf("Number of channels: %d\n",c)
% number of spectral windows
s = length(unique(data_id));
fprintf("Number of spectral windows: %d\n",s)
fprintf("Freqs:\n")
Freqs
fprintf("Ref freqs:\n")
RefFreq

sources = unique(field_id);
cyg = sources(end);

vis = vis(field_id==cyg,:);
flaging = flaging(field_id==cyg,:);
flag_row = flag_row(:,field_id==cyg)';
uvw = uvw(field_id==cyg,:);
weights_ch = weights_ch(:,field_id==cyg)';
times = times(:,field_id==cyg)';
data_id = data_id(:,field_id==cyg)';

bands = [8:240];

conf = 1;
nb_total = 0;
for i = 1 : s

    i

    Data{i} = vis(data_id==i-1,:);
    Flaging{i} = flaging(data_id==i-1,:);
    Flag_row{i} = flag_row(data_id==i-1,:);
    UVW{i} = uvw(data_id==i-1,:);
    Weights_ch{i} = weights_ch(data_id==i-1,:);
    Times{i} = times(data_id==i-1,:);

    for k = 1:length(bands)
        
        k
        
        j = bands(k);
        
        temp = Flaging{i}(:,j) | Flag_row{i}(:,:) | (Data{i}(:,j) == 0);

        ind = (temp==0);
        nb_total = nb_total + sum(ind);

        time_{k}{conf} = double(Times{i}(ind>0));

        I = Data{i}(:,j);
        w_ch = Weights_ch{i};
        nW{k}{conf} = double(w_ch(ind>0));
        y{k}{conf} = double(I(ind>0)) .* nW{k}{conf};
        f{k} = Freqs(i,j);

        u = UVW{i}(:,1);
        u = u(ind>0);

        v = - UVW{i}(:,2);
        v = v(ind>0);

        wave = speed/f{k};
        u = u ./ wave;
        v = v ./ wave;
        size(u)
        tmp = sqrt(u.^2 + v.^2);
        
        bmax_k(k,conf) = max(tmp(:));
        fprintf("ch:%d,bmax:%f\n",j,bmax_k(k,conf))

        u_{k}{conf} = u;
        v_{k}{conf} = v;

        pos{k}(conf) = length(y{k}{conf});
    end

end
fprintf("Number of total data points:%d\n",nb_total)

%% Setting imaging params
bmax = max(bmax_k(:));
bmin = min(bmax_k(:));
fprintf("Max baseline: %f\n", bmax)
fprintf("Min baseline: %f\n", bmin)
theta = 1 / (2*bmax);
theta_deg = theta * 180 / pi;
theta_arcsec = theta_deg * 3600;
fprintf("Theta:%f\n", theta_arcsec)

theta1 = 1 / (2*bmin);
theta_deg1 = theta1 * 180 / pi;
theta_arcsec1 = theta_deg1 * 3600;
fprintf("Theta:%f\n", theta_arcsec1)

pixel_size = 0.8;
dl = theta_arcsec / pixel_size;
fprintf("dl:%f\n", dl)

dl1 = theta_arcsec1 / pixel_size;
fprintf("dl:%f\n", dl1)

%% scale the uv-coverages w.r.t the desired dl
for i = 1 : length(u_)
%     fprintf("channel:%d\n",i)
    v_{i}{conf} = v_{i}{conf} * pi/(bmax*dl);
    u_{i}{conf} = u_{i}{conf} * pi/(bmax*dl);
%     fprintf("max |v|:%f\n",max(abs(v_{i}{conf})))
%     fprintf("max |u|:%f\n",max(abs(u_{i}{conf})))
end

nb_effch = floor(length(bands)/nb_aggch);
fprintf("Number of effective channels: %d\n", nb_effch)

for i = 1 : nb_effch
    fprintf("Effective channel:%d\n",i)
    u_temp = [];
    v_temp = [];
    time_temp = [];
    nW_temp = [];
    y_temp = [];
    for j = 1 : nb_aggch
        ind1 = nb_aggch*(i-1)+j;
        u_temp = [u_temp; u_{ind1}{conf}];
        v_temp = [v_temp; v_{ind1}{conf}];
        time_temp = [time_temp; time_{ind1}{conf}];
        nW_temp = [nW_temp; nW{ind1}{conf}];
        y_temp = [y_temp; y{ind1}{conf}];
    end
    if i == nb_effch
        rest_ch = length(bands) - nb_effch*nb_aggch;
        fprintf("Number of remaining channels:%d\n", rest_ch)
        for j = 1 : rest_ch
            ind1 = nb_aggch*nb_effch+j;
            u_temp = [u_temp; u_{ind1}{conf}];
            v_temp = [v_temp; v_{ind1}{conf}];
            time_temp = [time_temp; time_{ind1}{conf}];
            nW_temp = [nW_temp; nW{ind1}{conf}];
            y_temp = [y_temp; y{ind1}{conf}];
        end
        
    end
    u1{1}{1} = u_temp;
    v1{1}{1} = v_temp;
    time1{1}{1} = time_temp;
    nW1{1}{1} = nW_temp;
    y1{1}{1} = y_temp;
    save(['data_CGC_full_ch', num2str(i), '.mat'],'-v7.3','u1','v1','time1','nW1','y1') 
end

% save('data_CGC_full.mat','-v7.3','u_final','v_final','time_final','nW_final','y_final')


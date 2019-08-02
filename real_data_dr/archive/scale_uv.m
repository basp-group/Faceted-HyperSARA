function [u1, v1, pixel_size] = scale_uv(visibility_file_name,configuration_file_name, u1, v1)


speed = 299792458;

for conf = 1 : length(configuration_file_name) %visibility file
    conf
    k = 1;
    for m = 1 : length(visibility_file_name)  %visibility file
        
        m
        load(['CYG-' configuration_file_name{conf} '-' visibility_file_name{m} '-FULL.mat']);
        
        %% Split data
        % number of overall channels
        c = size(vis,2)
        % number of spectral windows
        s = length(unique(data_id));
        
        sources = unique(field_id);
        cyg = sources(end);
        
        uvw = uvw(field_id==cyg,:);
        flaging = flaging(field_id==cyg,:);
        flag_row = flag_row(:,field_id==cyg)';
        data_id = data_id(:,field_id==cyg)';
        
        for i = 1 : s
            
            i
            
            UVW{i} = uvw(data_id==i-1,:);
            Flaging{i} = flaging(data_id==i-1,:);
            Flag_row{i} = flag_row(data_id==i-1,:);
            
            for j = 1 : c
                
                temp = Flaging{i}(:,j) | Flag_row{i}(:,:);
                ind = (temp==0);
                
                f{k}{conf} = Freqs(i,j);
                
                u = UVW{i}(:,1);
                u = u(ind>0);
                
                v = - UVW{i}(:,2);
                v = v(ind>0);
                
                wave = speed/f{k}{conf};
                u = u ./ wave;
                v = v ./ wave;
                bmax_k(k,conf) = max(sqrt(u.^2 + v.^2));
                
                u_{k}{conf} = u;
                v_{k}{conf} = v;
                
                k = k + 1;
            end
            
        end
        
    end
end

%% Setting imaging params
bmax = max(bmax_k(:));
theta = 1 / (2*bmax);
theta_deg = theta * 180 / pi;
theat_arcsec = theta_deg * 3600

%dl = 2.5;
%pixel_size = round(theat_arcsec/dl,4)

pixel_size = 0.06  %we write it explicitly to be exact and copy the same value for wsclean
dl = theat_arcsec / pixel_size

%% scale the uv-coverages w.r.t the desired dl
for i = 1 : length(u1)
        v1{i} = v1{i} * pi/(bmax * dl);
        u1{i} = u1{i} * pi/(bmax * dl);
end












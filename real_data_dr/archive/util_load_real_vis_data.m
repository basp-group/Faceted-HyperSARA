function [y_final, u_final, v_final, nW_final, f, time_final, pos] = util_load_real_vis_data(visibility_file_name,configuration_file_name, ch_ind)


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
        
        vis = vis(field_id==cyg,:);
        flaging = flaging(field_id==cyg,:);
        flag_row = flag_row(:,field_id==cyg)';
        uvw = uvw(field_id==cyg,:);
        weights_ch = weights_ch(field_id==cyg,:);
        times = times(:,field_id==cyg)';
        data_id = data_id(:,field_id==cyg)';
        
        for i = 1 : s
            
            i
            
            Data{i} = vis(data_id==i-1,:);
            Flaging{i} = flaging(data_id==i-1,:);
            Flag_row{i} = flag_row(data_id==i-1,:);
            UVW{i} = uvw(data_id==i-1,:);
            Weights_ch{i} = weights_ch(data_id==i-1,:,:);
            Times{i} = times(data_id==i-1,:);
            
            for j = ch_ind
                
                temp = Flaging{i}(:,j) | Flag_row{i}(:,:);
                ind = (temp==0);
                
                time_{k}{conf} = double(Times{i}(ind>0));
                
                I = Data{i}(:,j);
                w_ch = Weights_ch{i}(:,j);
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
                bmax_k(k,conf) = max(sqrt(u.^2 + v.^2));
                
                u_{k}{conf} = u;
                v_{k}{conf} = v;
                
                pos{k}(conf) = length(y{k}{conf});
                
                k = k + 1;
            end
            
        end
        
    end
end

for i = 1 : length(u_)
    u_temp = [];
    v_temp = [];
    time_temp = [];
    nW_temp = [];
    y_temp = [];
    for j = 1 : length(u_{i}) % aggregate configurations once done (blocking per config): simple reshape
        u_temp = [u_temp; u_{i}{j}];
        v_temp = [v_temp; v_{i}{j}];
        time_temp = [time_temp; time_{i}{j}];
        nW_temp = [nW_temp; nW{i}{j}];
        y_temp = [y_temp; y{i}{j}];
    end
    u_final{i} = u_temp;
    v_final{i} = v_temp;
    time_final{i} = time_temp;
    nW_final{i} = nW_temp;
    y_final{i} = y_temp;
end


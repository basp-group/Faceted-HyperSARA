function [y_final, u_final, v_final, nW_final, f_final, time_final, pos_final] = util_load_real_vis_data_dr(visibility_file_name, nSpw, ind)

% nSpw: total number of sepctral windows (i.e., over the total number of 
% configurations considered)
% pos: keep track of the different configurations
% 32 problems in total: 16 spw per MS for each config, 2 MS files for each
% configuration

% see if anything else is needed at this stage... (one point unclear from 
% Abdullah's explanations)

y_final = [];
u_final = [];
v_final = [];
nW_final = [];
f_final = [];
time_final = [];
pos_final = [];

for n = 1:nSpw
    filename = [visibility_file_name, num2str(n), '.mat'];
    file = matfile(filename);
    f_final = [f_final; cell2mat(file.f(1,ind))];
    nW_final = [nW_final; cell2mat(file.nWw(1,ind))];
    pos_final = [pos_final; cell2mat(file.pos(1,ind))];
    time_final = [time_final; cell2mat(file.time_(1,ind))];
    u_final = [u_final; cell2mat(file.uw(1,ind))];
    v_final = [v_final; cell2mat(file.vw(1,ind))];    
    y_final = [y_final; cell2mat(file.y(1,ind))];
end

end

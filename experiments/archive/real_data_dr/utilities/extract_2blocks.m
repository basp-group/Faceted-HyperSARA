% from 9 data blocks to 2 data blocks (a, c configuration)

nChannels = 32;
nBlocks = 2;
c_config_block_index = 2;

new_file_y = matfile('/home/basphw/mjiang/Data/pthouvenin/basp_sharing/Ming/data_mnras_dr/CYG_y.mat');
new_file_u = matfile('/home/basphw/mjiang/Data/pthouvenin/basp_sharing/Ming/data_mnras_dr/CYG_u.mat');
new_file_v = matfile('/home/basphw/mjiang/Data/pthouvenin/basp_sharing/Ming/data_mnras_dr/CYG_v.mat');
new_file_nW = matfile('/home/basphw/mjiang/Data/pthouvenin/basp_sharing/Ming/data_mnras_dr/CYG_nW.mat');

new_file_2b_y = matfile('/home/basphw/mjiang/Data/mjiang/extract_real_data//CYG_2b_y.mat', 'Writable', true);
new_file_2b_u = matfile('/home/basphw/mjiang/Data/mjiang/extract_real_data/CYG_2b_u.mat', 'Writable', true);
new_file_2b_v = matfile('/home/basphw/mjiang/Data/mjiang/extract_real_data/CYG_2b_v.mat', 'Writable', true);
new_file_2b_nW = matfile('/home/basphw/mjiang/Data/mjiang/extract_real_data/CYG_2b_nW.mat', 'Writable', true);

for l = 1:nChannels
    y_tmp = new_file_y.y(l, 1);
    u_tmp = new_file_u.u(l, 1);
    v_tmp = new_file_v.v(l, 1);
    nW_tmp = new_file_nW.nW(l, 1);

    y_tmp1 = cell(nBlocks, 1);
    u_tmp1 = cell(nBlocks, 1);
    v_tmp1 = cell(nBlocks, 1);
    nW_tmp1 = cell(nBlocks, 1);

    oldBlocks = length(y_tmp{1});

    b_pos = 1;
    flag = 1;
    for m = 1:nBlocks
        for n = b_pos:oldBlocks
            y_tmp1{m, 1} = [y_tmp1{m, 1}; y_tmp{1}{n}];
            u_tmp1{m, 1} = [u_tmp1{m, 1}; u_tmp{1}{n}];
            v_tmp1{m, 1} = [v_tmp1{m, 1}; v_tmp{1}{n}];
            nW_tmp1{m, 1} = [nW_tmp1{m, 1}; nW_tmp{1}{n}];
            b_pos = n + 1;
            if b_pos > c_config_block_index && flag
                flag = 0;
                break
            end
        end
    end
    new_file_2b_y.y(l, 1) = {y_tmp1};
    new_file_2b_u.u(l, 1) = {u_tmp1};
    new_file_2b_v.v(l, 1) = {v_tmp1};
    new_file_2b_nW.nW(l, 1) = {nW_tmp1};
end

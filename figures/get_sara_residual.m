function get_sara_residual(pathname, src_filename, dest_filename, Ny, Nx, nChannels)

        res = zeros(Ny, Nx, nChannels);
        for l = 1:nChannels
            m = matfile(strcat(src_filename(l), ".mat"));
            res(:,:,l) = m.res; % see if ok
        end
        fitswrite(res, fullfile(pathname, strcat(dest_filename, ".fits")));
end

% get_sara_residual(@(l) strcat("../results/test/test_cygASband_Cube_512_1024_20_sara_none_srf=2_Ny=512_Nx=1024_L=20_Qy=1_Qx=1_Qc=20_ind=",num2str(l),"_gam=1_gambar=1_overlap=0_0_rw_type=dirty_snr=40.mat"), "test_sara", 512, 1024, 20)
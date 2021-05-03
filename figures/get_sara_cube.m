function get_sara_cube(pathname, Ny, Nx, nChannels, gam, homotopy, rwtype)

    filename = @(channel) strcat(...
        "x_cygASband_Cube_512_1024_20_sara_srf=2_none_Qy=1_Qx=1_Qc=20_ind=", ...
        num2str(channel), "gam=", num2str(gam), ...
        "_homotopy=1_rwtype=", rwtype,"homotopy=", homotopy, "_snr=40");

    cubename = @(channel) strcat(...
        "x_cygASband_Cube_512_1024_20_sara_srf=2_none_Qy=1_Qx=1_Qc=20_", ...
        "gam=", num2str(gam), "_homotopy=1_rwtype=", rwtype, ...
        "homotopy=", homotopy, "_snr=40");

    % x_cygASband_Cube_512_1024_20_sara_srf=2_none_Qy=1_Qx=1_Qc=20_ind=7_gam=10_homotopy=1_rwtype=dirty_snr=40.fits

    x = zeros(Ny, Nx, nChannels);

    for l = 2:nChannels
        fullFilename = fullfile(pathname, strcat(filename(l), ".fits"));
        x_ = fitsread(fullFilename);
        x(:,:,l) = x_;
    end
    fitswrite(x, cubename);
end
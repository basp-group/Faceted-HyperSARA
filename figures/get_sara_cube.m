function get_sara_cube(pathname, Ny, Nx, nChannels, gam, homotopy, rwtype, ...
    pathToGrountruth)

    filename = @(channel) strcat(...
        "x_cygASband_Cube_512_1024_20_sara_srf=2_none_Qy=1_Qx=1_Qc=20", ...
        "_ind=", num2str(channel), "_gam=", num2str(gam), ...
        "_homotopy=", num2str(homotopy), "_rwtype=", rwtype, "_snr=40");

    cubename = strcat("x_cygASband_Cube_512_1024_20_sara_srf=2_none_Qy=1_Qx=1_Qc=20_gam=", ... 
    num2str(gam), "_homotopy=", num2str(homotopy), "_rwtype=", rwtype, ...
    "_snr=40.fits");

    % x_cygASband_Cube_512_1024_20_sara_srf=2_none_Qy=1_Qx=1_Qc=20_ind=7_gam=10_homotopy=1_rwtype=dirty_snr=40.fits

    x = zeros(Ny, Nx, nChannels);

    for l = 1:nChannels
        disp(fullfile(pathname, strcat(filename(l), '.fits')));
        fullFilename = fullfile(pathname, strcat(filename(l), '.fits'));
        x_ = fitsread(fullFilename);
        x(:,:,l) = x_;
    end
    disp(fullfile(pathname,cubename))
    fitswrite(x, fullfile(pathname,cubename));

    % evaluate full and per channel reconstruction snr
    x0 = fitsread(pathToGrountruth);
    err_per_channel = squeeze(sum(x0.^2, [1,2])./sum((x - x0).^2, [1,2]));
    SNR = 10*log10(sum(err_per_channel));
    SNR_per_channel = 10*log10(err_per_channel);

    fprintf("SNR = %e \n", SNR);
    fprintf("SNR per channel = %e \n\n", SNR_per_channel);

end
clc; clear all; close all;
format compact

% algo = ["hypersara", "cw"];
% gam = ["1e-1","1","10"];
% gambar = ["1e-1","1","10"];

algo = ["hypersara", "cw"];
gam = [1e-1, 1, 10];
gam_bar = [1e-1, 1, 10];

src_filename = @(gam, gambar, algo) strcat("test_cygASband_Cube_512_1024_20_", algo, "_triangular_srf=2_Ny=512_Nx=1024_L=20_Qy=2_Qx=2_Qc=1_ind=1_gam=", num2str(gam), "_gambar=", num2str(gambar), "_overlap=0.5_0.5_rw_type=dirty_snr=40.mat");

resname = @(gam, gambar, algo) strcat("res_", algo, "_triangular_srf=2_Ny=512_Nx=1024_L=20_Qy=2_Qx=2_Qc=1_ind=1_gam=", num2str(gam), "_gambar=", num2str(gambar), "_overlap=0.5_0.5_rw_type=dirty_snr=40.mat");

imname = @(gam, gambar, algo) strcat("x_", algo, "_triangular_srf=2_Ny=512_Nx=1024_L=20_Qy=2_Qx=2_Qc=1_ind=1_gam=", num2str(gam), "_gambar=", num2str(gambar), "_overlap=0.5_0.5_rw_type=dirty_snr=40.mat");

for k = 1:numel(gam)
    for l = 1:numel(gam_bar)
        for m = 1:numel(algo)
            f = matfile(src_filename(gam(k), gambar(l),algo(m)));
            fitswrite(f.res, strcat(resname(gam(k), gambar(l),algo(m)), ".fits"));
            fitswrite(f.xsol, strcat(imname(gam(k), gambar(l),algo(m)), ".fits"));
        end
    end
end

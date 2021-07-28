flag_generateUndersampledCube = false;
flag_generateCoverage = false;

%% Generate/load ground-truth image cube
% reference_cube_path = strcat(data_path, image_name, '.fits');
% if flag_generateCube
%     % frequency bandwidth from 1 to 2 GHz
%     f = linspace(1,2,nChannels);
%     emission_lines = 0; % insert emission lines on top of the continuous spectra
%     % [x0,X0] = Generate_cube(reference_cube_path,f,emission_lines);
%     [x0,X0] = Generate_cube_W28(reference_cube_path,f,emission_lines);
%     [Ny, Nx, nChannels] = size(x0);
%     if flag_generateUndersampledCube
%         % undersampling factor for the channels
%         unds = 4; % take 1/unds images
%         [x0,X0,f,nChannels] = Generate_undersampled_cube(x0,f,Ny,Nx,nChannels,unds);
%     end
%     fitswrite(X0.', strcat(cube_path, '.fits'));
%     fitsdisp(strcat(cube_path, '.fits'));
% else
%     X0 = fitsread(strcat(cube_path, '.fits')).';
%     Nx = sqrt(size(X0, 1));
%     Ny = Nx;
%     nChannels = size(X0, 2);
%     x0 = reshape(X0, [Ny, Nx, nChannels]);
% end
% % frequency bandwidth from 1 to 2 GHz
% f = linspace(1, 2, nChannels);


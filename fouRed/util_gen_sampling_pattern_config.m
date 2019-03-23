%% sampling pattern parameters
% options 'gaussian', 'file', 'gaussian+large-holes', 'file+undersample', ''gaussian+missing-pixels'
% sampling_pattern = 'file+undersample';
% sampling_pattern = 'gaussian';
sampling_pattern = 'gaussian+large-holes';
% sampling_pattern = 'gaussian+missing-pixels';
 
% param_sampling.file_name = 'data/meerkat2h.ar.uvw.dat'; % file name for uv coverage
% param_sampling.file_name = '~/Data/MeasSets/ska.2h.ar.uvw.dat'; % file name for uv coverage
% param_sampling.N=256*256;
% param_sampling.p = 50; % number of measurements as proportion of number of pixels to recover
param_sampling.p = ceil(visibSize/N);
param_sampling.hole_number = 4000; % number of holes to introduce for 'gaussian+large-holes'
param_sampling.hole_prob = 0.1; % probability of single pixel hole for 'gaussian+missing-pixels'
param_sampling.hole_size = pi/60; % size of the missing frequency data
param_sampling.fpartition = [pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% param_sampling.fpartition = [0, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% param_sampling.fpartition = [-0.25*pi, 0, 0.25*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% param_sampling.fpartition = [-0.3*pi, -0.15*pi, 0, 0.15*pi, 0.3*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% param_sampling.fpartition = [-0.35*pi, -0.25*pi, -0.15*pi, 0, 0.15*pi, 0.25*pi, 0.35*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
param_sampling.sigma = pi/4; % variance of the gaussion over continous frequency
param_sampling.sigma_holes = pi/3; % variance of the gaussion for the holes

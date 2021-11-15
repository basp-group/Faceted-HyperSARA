clear all; close all; format compact
%%
% https://fr.mathworks.com/matlabcentral/answers/474980-extract-info-from-json-file-by-matlab

fileName = 'setup_test.json'; % filename in JSON extension
fid = fopen(fileName); % Opening the file
raw = fread(fid,inf); % Reading the contents
str = char(raw'); % Transformation
fclose(fid); % Closing the file
data = jsondecode(str); % Using the jsondecode function to parse JSON from string
param_nufft = data{2, 1}.nufft;

% testing parallel random generation
seed = 1;
n = 10;
[stream{1:n}] = RandStream.create('threefry4x64_20','NumStreams',n,'Seed',seed);
parfor k = 1:n
    r = randn(stream{k},[1 3]);
    disp(r);
end

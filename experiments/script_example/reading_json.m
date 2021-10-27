clear all; close all; format compact
%%
% https://fr.mathworks.com/matlabcentral/answers/474980-extract-info-from-json-file-by-matlab

fileName = 'setup.json'; % filename in JSON extension
fid = fopen(fileName); % Opening the file
raw = fread(fid,inf); % Reading the contents
str = char(raw'); % Transformation
fclose(fid); % Closing the file
data = jsondecode(str); % Using the jsondecode function to parse JSON from string
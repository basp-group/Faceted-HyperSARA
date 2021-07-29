function [v2, norm_res, t_block, proj] = initiliaze_data_worker(y)

nChannels = numel(y);
norm_res = cell(nChannels, 1);
v2 = cell(nChannels, 1);
t_block = cell(nChannels, 1);
proj = cell(nChannels, 1);

% loop over channels
for l = 1:nChannels
    norm_res{l} = cell(numel(y{l}),1);
    v2{l} = cell(numel(y{l}),1);
    t_block{l} = cell(numel(y{l}),1);
    proj{l} = cell(numel(y{l}),1);
    % loop over data blocks
    for b = 1 : numel(y{l})
        norm_res{l}{b} = norm(y{l}{b});
        v2{l}{b} = zeros(numel(y{l}{b}) ,1);
        t_block{l}{b} = 0;
        proj{l}{b} = zeros(numel(y{l}{b}), 1);
    end
end

end

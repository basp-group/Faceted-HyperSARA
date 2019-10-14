nSpw = 16;          % number of spectral channels per MS file
nChannels = 2*nSpw; % total number of "virtual" channels (i.e., after
% concatenation) for the real dataset considered
nBlocks = 2;        % number of data blocks (needs to be known beforehand,
% quite restrictive here), change l.70 accordingly

for i = 1:nChannels
    i
    tmp = load(['/lustre/home/shared/sc004/dr_2b_result_real_data/CYG_DR=', num2str(i), '.mat'], 'H');
    H = tmp.H{1};
    for j = 1:nBlocks
        a = H{j};
        H{j} = a(:);
        disp("Block " + j + " vectorised");
    end
    save(['/lustre/home/shared/sc004/dr_2b_result_real_data/CYG_DR=', num2str(i), '.mat'], 'H', '-append');
    disp("Channel " + i + " finished");
end
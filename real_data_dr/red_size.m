function red_size(ch, realdatablocks, reduction_version, fouRed_gamma)

% diaryFname = ['diary_', num2str(ch(1)), '_', num2str(ch(end)), '_', num2str(realdatablocks), 'b_fouRed',...
%         num2str(reduction_version), '_perc', num2str(fouRed_gamma), '.txt'];
% if exist(diaryFname, 'file')
%     delete(diaryFname)
% end
% diary(diaryFname)

fprintf('Channel number %d\n', ch);
fprintf('Reduction version %d\n', reduction_version);
fprintf('Data blocks: %d\n', realdatablocks);
fprintf('Reduction level: %f sigma\n', fouRed_gamma);

% ch = 1:32;
% realdatablocks = 2;
% reduction_version = 2;
% fouRed_gamma = 20;
nChannels = length(ch);

raw_ctotal = 0;
raw_atotal = 0;
raw_total = 0;
nb_ctotal = 0;
nb_atotal = 0;
nb_total = 0;
for i = 1:nChannels 
    ch(i)
    tmp = load(['/lustre/home/shared/sc004/dr_', num2str(realdatablocks), 'b_result_real_data/CYG_DR_cal_', num2str(realdatablocks), 'b_ind6_fouRed',...
        num2str(reduction_version), '_th', num2str(fouRed_gamma),'=', num2str(ch(i)), '.mat'], 'H', 'Wm');
%     H{i,1} = tmp.H{1,1};
%     W{i,1} = tmp.W{1,1};
%     yT{i,1} = tmp.yT{1,1};
%     T{i,1} = tmp.T{1,1};
%     Wm{i,1} = tmp.Wm{1,1};
    nb = 0;
    raw = 0;
    for j = 1:realdatablocks
        nb_b = sum(tmp.Wm{1,1}{j});
        raw_b = numel(tmp.Wm{1,1}{j});
        fprintf('Number of points of block %d of channel %d: %d\n', j, i, nb_b)
        if (j == 1 && realdatablocks == 2) || ((j == 1||j == 2) && realdatablocks == 9)
            nb_ctotal = nb_ctotal + nb_b;
            raw_ctotal = raw_ctotal + raw_b;
        elseif (j == 2 && realdatablocks == 2) || (j > 2 && realdatablocks == 9)
            nb_atotal = nb_atotal + nb_b;
            raw_atotal = raw_atotal + raw_b;
        end
        nb = nb + nb_b;
        raw = raw + raw_b;
    end
    raw_total = raw_total + raw;
    nb_total = nb_total + nb;
    fprintf('Number of raw data points of channel %d: %d\n', i, raw)
    fprintf('Number of data points of channel %d: %d\n', i, nb)
    fprintf('H memory of channel %d\n', i)
    a = tmp.H{1,1};
    whos a
    clear tmp a
end
fprintf('Total number of raw C-block data points: %d\n', raw_ctotal)
fprintf('Total number of raw A-block data points: %d\n', raw_atotal)
fprintf('Total number of raw data points: %d\n', raw_total)

fprintf('Total number of C-block data points: %d\n', nb_ctotal)
fprintf('Total number of A-block data points: %d\n', nb_atotal)
fprintf('Total number of data points: %d\n', nb_total)
% diary off
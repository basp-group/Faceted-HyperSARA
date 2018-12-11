function [y, Sigmas, masks] = util_gen_sing_block_structure(x, Sigma, mask, param)

%% optional input arguments
if ~isfield(param, 'use_uniform_partitioning'), param.use_uniform_partitioning = 1; end
if ~isfield(param, 'uniform_partitioning_no'), param.uniform_partitioning_no = 1; end

if param.use_uniform_partitioning
    nb_sing = length(Sigma);
    step = floor(nb_sing/param.uniform_partitioning_no);
    masks = cell(param.uniform_partitioning_no,1);
    [masks{:}] = deal(sparse(false(length(mask),1)));        % Initialization, each cell element is a vector of image size
    beg_ind = 1;                                     % Beginning index of the mask under processing
    for i = 1:param.uniform_partitioning_no-1
        tmp = find(mask(beg_ind:end), step, 'first');    % Find non-zero elements of the mask which have not been treated, the number of non-zero elements corresponds to the size of block.
        end_ind = tmp(end);                          % Ending non-zero index of the sub-mask
        masks{i,1}(beg_ind:beg_ind+end_ind-1) = mask(beg_ind:beg_ind+end_ind-1);
        beg_ind = beg_ind + end_ind;                      % Update the beginning index
        ind = (i-1)*step+1:i*step;
        Sigmas{i} = Sigma(ind);
        y{i} = x(ind);
    end
    masks{param.uniform_partitioning_no,1}(beg_ind:end) = mask(beg_ind:end);
    ind = (param.uniform_partitioning_no-1)*step+1:nb_sing;
    Sigmas{param.uniform_partitioning_no} = Sigma(ind);
    y{param.uniform_partitioning_no} = x(ind);
end
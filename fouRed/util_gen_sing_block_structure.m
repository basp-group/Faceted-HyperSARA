function [y, nss, Sigmas, masks] = util_gen_sing_block_structure(x, ns, Sigma, mask, param)

%% optional input arguments
if ~isfield(param, 'use_uniform_partitioning'), param.use_uniform_partitioning = 1; end
if ~isfield(param, 'uniform_partitioning_no'), param.uniform_partitioning_no = 1; end

if param.use_uniform_partitioning
    nb_sing = length(Sigma);
    step = floor(nb_sing/param.uniform_partitioning_no);    
    beg_ind = 1;                                     % Beginning index
    mask_ind = find(mask);                           % Index vector of singular values
    for i = 1:param.uniform_partitioning_no-1
        masks{i} = mask_ind(beg_ind:beg_ind+step-1);
        Sigmas{i} = Sigma(beg_ind:beg_ind+step-1);
        y{i} = x(beg_ind:beg_ind+step-1);
        nss{i} = ns(beg_ind:beg_ind+step-1);
        beg_ind = beg_ind + step;
    end
    % The last partition may contain different number of elements
    masks{param.uniform_partitioning_no} = mask_ind(beg_ind:end);
    Sigmas{param.uniform_partitioning_no} = Sigma(beg_ind:end);
    y{param.uniform_partitioning_no} = x(beg_ind:end);
    nss{param.uniform_partitioning_no} = ns(beg_ind:end);
    
elseif param.use_sort_uniform_partitioning
    nb_sing = length(Sigma);
    step = floor(nb_sing/param.uniform_partitioning_no);  
    beg_ind = 1;                                     % Beginning index
    mask_ind = find(mask);                           % Index vector of singular values
    [Sigma_sort, ind_sort] = sort(abs(Sigma));
    mask_ind_sort = mask_ind(ind_sort);
    x_sort = x(ind_sort);
    ns_sort = ns(ind_sort);
    for i = 1:param.uniform_partitioning_no-1
        masks{i} = mask_ind_sort(beg_ind:beg_ind+step-1);
        Sigmas{i} = Sigma_sort(beg_ind:beg_ind+step-1);
        y{i} = x_sort(beg_ind:beg_ind+step-1);
        nss{i} = ns_sort(beg_ind:beg_ind+step-1);
        beg_ind = beg_ind + step;
    end
    % The last partition may contain different number of elements
    masks{param.uniform_partitioning_no} = mask_ind_sort(beg_ind:end);
    Sigmas{param.uniform_partitioning_no} = Sigma_sort(beg_ind:end);
    y{param.uniform_partitioning_no} = x_sort(beg_ind:end);
    nss{param.uniform_partitioning_no} = ns_sort(beg_ind:end);
end
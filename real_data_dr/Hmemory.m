precision = 1e-16;

sigma = 1;
p = normcdf([-sigma sigma]);
prob = p(2) - p(1);

HRed = cell(2,1);

for j = 1:2
    peak = max(max(abs(H{1}{j})));
    HRed{j} = H{1}{j} .* (abs(H{1}{j}) > peak * precision);
    d_mat = full(abs(diag(HRed{j})));
    
%     Mask = (d_mat > 0);
    d_mat_sort = sort(d_mat);
    d_mat_sort_cumsum = cumsum(d_mat_sort);
    d_mat_sum = d_mat_sort_cumsum(end); % whole energy of d_mat
    th_energy = d_mat_sum * (1 - prob); % threshold according to the k-sigma rule
    th = d_mat_sort(find(d_mat_sort_cumsum >= th_energy, 1, 'first'));
    Mask = (d_mat >= th);
    
    HRed{j} = HRed{j}(Mask,:);
end

clear HRed
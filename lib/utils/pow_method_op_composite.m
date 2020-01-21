function val = pow_method_op_composite(Hp, Wp, Ap, Atp, Tp, aWp, Q, K, nchannel_per_worker, im_size)
%Computes the maximum eigen value of the compund 
%operator AtA
%
% > A       composite lambda function (to be defined as a constant on each worker)
% > At      composite lambda function for the full operator (per channel)
% > Q       number of facets (workers)
% > K       number of data nodes
% > nchannel_per_worker  number of channels per data worker [K, 1]
% > im_size image size (spatial dimensions) [2, 1]
%
% P.-A. Thouvenin, [04/01/2020]

p = 1 + 10^(-6);
pnew = 1;
epsilon = 10^(-8);
n = 1;
nmax = 200;
cond = abs(pnew - p) / pnew;

spmd
    if labindex > Q
         % each worker owns a portion of the channels
        x = randn([im_size, nchannel_per_worker(labindex-Q)]);
        norm2_x = norm(x(:))^2;
    end
end
norm_x = 0;
for i=Q+1:Q+K
    norm_x = norm_x + norm2_x{i};
end
norm_x = sqrt(norm_x);

% Iterations
while (cond >= epsilon && n < nmax)
    spmd
        if labindex > Q 
%            xnew = At(A(x/norm_x)); % input normalised
           xtmp = HS_operatorGtPhi(x/norm_x, Hp, Wp, Ap, Tp, aWp);
           xnew = HS_operatorGtPhi_t(xtmp, Hp, Wp, Atp, Tp, aWp);
                    
           squared_norm2_diff = norm(xnew(:) - x(:))^2;
           
           x = xnew;
           norm2_x = norm(x(:))^2;
        end
    end
    
    norm_x = 0;
    squared_norm_diff = 0;
    for i=Q+1:Q+K
        norm_x = norm_x + norm2_x{i};
        squared_norm_diff = squared_norm_diff + squared_norm2_diff{i};
    end
    norm_x = sqrt(norm_x);
    squared_norm_diff = sqrt(squared_norm_diff);
    
    cond = squared_norm_diff / norm_x;
    n = n+1 ;   
    display(['norm it n=' num2str(n) ' norm=' num2str(norm_x)])
end

val = norm_x; % output not normalised

end


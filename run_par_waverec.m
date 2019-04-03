function [v1_, u1_, l21_] = run_par_waverec(v1_, x_, weights1_, beta1)

l21_ = 0;
for i = 1 : size(x_, 3)
    r1 = v1_(:, i) +  Psit(x_(:, :, i));
    
    
    
    l2 = sqrt(sum(abs(r1).^2,2));
    l2_soft = max(l2 - beta1*weights1_{i}, 0)./(l2+eps);
    v1_{i} = r1 - l2_soft.*r1;
    u1_{i} = Psi(v1_{i});
    
    % local L21 norm of current solution
    l21_ = l21_ + norm(l2(:),1);
end
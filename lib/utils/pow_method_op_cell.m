function val = pow_method_op_cell(A, At, im_size)
%Computes the maximum eigen value of the compund 
%operator AtA
% 
N = im_size(1);
c = im_size(2); 

for i = 1 : c
    x{i}=randn(N);
end

norm_x = norm(cell2mat(x'));
for i = 1 : c
    x{i}=x{i}/norm_x;
end

p = 1 + 10^(-6) ;
pnew = 1 ;

n = 1 ;

epsilon = 10^(-8) ;

nmax = 200;

cond = abs( pnew-p ) / pnew ;

% Iterations

while ( cond >= epsilon && n < nmax)
    xnew=At(A(x));
    p=pnew;
    pnew=norm(cell2mat(xnew')) / norm(cell2mat(x')) ;
        
    cond = abs(  pnew-p ) / pnew ;
    
%     pp(n) = pnew ;
    
    x = xnew;

norm_x = norm(cell2mat(x'));
for i = 1 : c
    x{i}=x{i}/norm_x;
end

    n = n+1 ;
    
    display(['norm it n=' num2str(n) ' norm=' num2str(p)])
end
% figure
% plot(pp)
% pause

val = p ;
% x=randn(im_size);
% x=x/norm(x(:));
% init_val=1;
% 
% for k=1:max_iter
%     x=At(A(x));
%     val=norm(x(:));
%     rel_var=abs(val-init_val)/init_val;
%     if (verbose > 0)
%         fprintf('Iter = %i, norm = %e \n',k,val);
%     end
%     if (rel_var < tol)
%         break;
%     end
%     init_val=val;
%     x=x/val;
%     
% end


end


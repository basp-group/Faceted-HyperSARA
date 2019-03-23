function val = pow_method_op(A, At, im_size)
%Computes the maximum eigen value of the compund 
%operator AtA
%   
x=randn(im_size);
x=x/norm(x(:));

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
    pnew=norm(xnew(:)) / norm(x(:)) ;
    
    cond = abs(  pnew-p ) / pnew ;
    
%     pp(n) = pnew ;
    
    x = xnew;
    x=x/norm(x(:));

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


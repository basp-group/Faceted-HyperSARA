function val = pow_method_op(A, At, im_size)
%Computes the maximum eigen value of the compund 
%operator AtA
%   
x=randn(im_size);
x_norm=x/norm(x(:));

p = 1 + 10^(-6) ;
pnew = 1;

n = 1 ;

epsilon = 10^(-8) ;

nmax = 200;

cond = abs( pnew-p ) / pnew ;

% Iterations

while ( cond >= epsilon && n < nmax)
    
    xnew = At(A(x_norm)); % input normalised
    
    cond = norm(xnew(:) - x(:))  / norm(x(:));
        
    x = xnew;
    x_norm =x/norm(x(:));

    n = n+1 ;
    
    display(['norm it n=' num2str(n) ' norm=' num2str(norm(xnew(:)))])
end

val = norm(xnew(:)); % output not normalised

end


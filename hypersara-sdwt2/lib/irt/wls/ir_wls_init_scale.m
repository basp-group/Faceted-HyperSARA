 function xnew = ir_wls_init_scale(A, y, xold)
%function xnew = ir_wls_init_scale(A, y, xold)
%|
%| Find the scaled version of the image "x" that best fits the data.
%|
%| out
%|	xnew = scale * xold
%|	where scale = argmin_s |y - s A x|_2^2
%|
%| 2014, Jeff Fessler, University of Michigan

if nargin < 3
	xold = A' * y; % default that is sensible only for some applications
end
tmp = A * xold;
scale = sum(col(conj(tmp) .* y), 'double') / sum(col(abs(tmp).^2), 'double');
xnew = scale * xold;

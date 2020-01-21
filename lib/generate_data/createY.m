function Y = createY(y, N)
%createY Create an Hermitian matrix Y from the coefficients of its lower 
%   half.
%
% Input:
% > y  coefficients of the lower part of the Hermitian matrix [N*(N-1)/2,1]
% > N  size of the Hermitian matrix
%
% Output:
% < Y  Hermitian matrix composed of the entries given in y [N,N]
%-------------------------------------------------------------------------%
% P.-A. Thouvenin, October 11 2017.
%-------------------------------------------------------------------------%
%%
id = getIndicesLow(N); % indices of the lower part of Y
Y = zeros(N);
Y(id) = conj(y(:));
Y = Y + Y';            % build the Hermitian matrix

end

 function [y1, y2, y3] = solver3(a, b, c, d)
%function [y1, y2, y3] = solver3(a, b, c, d)
%
% The function solves the 3rd-degree equation:
%
%        a * x^3 + b * x^2 + c * x + d = 0
%
% and returns the 3 roots in the output vectors y1, y2, y3.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (27-04-2017)
% Author  : Giovanni Chierchia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2017
%
% This file is part of the codes provided at http://proximity-operator.net
%
% By downloading and/or using any of these files, you implicitly agree to 
% all the terms of the license CeCill-B (available online).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% check inputs
if any( a(:)==0 )
    error('A cubic equation is expected');
end
%-----%


% reduction to a depressed cubic (t^3 + p*t + q = 0)
p = c ./ a - b.^2 ./ (3*a.^2);
q = (2*b.^3/27 - a.*b.*c/3 + a.^2.*d) ./ a.^3;

% discriminant
DD = (p/3).^3 + (q/2).^2;

% output vectors
y1 = zeros( size(DD) );
y2 = zeros( size(DD) );
y3 = zeros( size(DD) );

% 1st case: 3 real unequal roots
idx = DD < 0;
phi = acos( -q(idx)/2 ./ sqrt(abs(p(idx).^3/27)) );
tau = 2 * sqrt( abs(p(idx)/3) );
y1(idx) =  tau .* cos( phi    /3);
y2(idx) = -tau .* cos((phi+pi)/3);
y3(idx) = -tau .* cos((phi-pi)/3);

%2nd case: 1 real root + 2 conjugate complex roots
idx = DD > 0;
z1 = -q(idx)/2;
z2 = sqrt( DD(idx) );
u = nthroot(z1+z2, 3);
v = nthroot(z1-z2, 3);
y1(idx) = u + v;
e1 = (-1 + 1i * sqrt(3)) / 2;
e2 = (-1 - 1i * sqrt(3)) / 2;
y2(idx) = u*e1 + v*e2;
y3(idx) = u*e2 + v*e1;

% 3rd case: 1 simple real root + 1 double real root
idx = DD == 0;
y1(idx) =  3.0 * q(idx) ./ p(idx);
y2(idx) = -1.5 * q(idx) ./ p(idx);
y3(idx) = y2(idx);

% correction to the original cubic
t = b ./ (3*a);
y1 = y1 - t;
y2 = y2 - t;
y3 = y3 - t;
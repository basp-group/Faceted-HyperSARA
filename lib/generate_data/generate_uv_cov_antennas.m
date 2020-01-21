function [u,v,w,na] = generate_uv_cov_antennas(antenna_position, x0, h, lat, dec, nstep)

% create u-v position of the  antennas
na = length(antenna_position); % number of antennas

% compute baseline lengths
b = zeros(na,2) ;
b(:,1) = sqrt( (antenna_position(:,1)-x0(1)).^2 + (antenna_position(:,2)-x0(2)).^2 );
b(:,2) = atan2( antenna_position(:,1)-x0(1), antenna_position(:,2)-x0(2) );

% selection in the u-v plane
u = zeros(nstep, na);
v = zeros(nstep, na);
w = zeros(nstep, na);
for a = 1:na
    X = uvTrack(h, b(a, 1), b(a, 2), 0., lat, dec, nstep);
    u(:, a) = X(:, 1);
    v(:, a) = X(:, 2);
    w(:, a) = X(:, 3);
end

end

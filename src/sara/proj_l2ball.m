function p = proj_l2ball(x, eps, y)

% projection of x onto the l2 ball centered in y with radius eps
p = x-y;
p = p* min(eps/norm(p(:)),1);
p = p+y;

end
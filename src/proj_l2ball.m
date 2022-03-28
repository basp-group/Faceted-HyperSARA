function p = proj_l2ball(x, eps, y)
% Compute the projection of the input vector :math:`x` onto the
% :math:`\ell_2` ball centered in :math:`y` with radius :math:`\varepsilon`.
%
% Parameters
% ----------
% x : array
%     Input vector.
% eps : double
%     Radius of the :math:`\ell_2` ball.
% y : array
%     Center of the :math:`\ell_2` ball.

% Returns
% -------
% p : array
%     Projection of :math:`x` onto the :math:`\ell_2` ball
%     :math:`\mathcal{B} (y, \varepsilon)`.

p = x - y;
p = p * min(eps / norm(p(:)), 1);
p = p + y;

end

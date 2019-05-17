function [w, wg] = generate_2D_gaussian(Nfacet, Npadded, sigx, sigy, bool)
% \Sigma = [sig(1), rho*prod(sig); rho*prod(sig), sig(2)]
% parameter: (10, 10, 0) (find angular parameterization fo the Gaussian distribution, use the same parameters as in ddfacet)

if ~bool
    N = Npadded;
    w = ones(Nfacet);
    % N = Npadded - Nfacet + 1;

    % Generate the 2D Gaussian window
    alpha_x = 2*(N(2) - 1)*sigx;
    alpha_y = 2*(N(1) - 1)*sigy;
    wy = gausswin(N(1), alpha_y);
    wx = gausswin(N(2), alpha_x); % possibly take a different window...
    wg = wy.*(wx.');

    % Convolve with a unit window of the facet size
    w = conv2(wg, w, 'same');

    % Threshold the smallest coefficients to zero / renormalize the window
    w(w < 1e-5) = 0; % set threshold as a parameter
    w = w/max(w(:)); % normalized window
else
    N = Npadded;
    w = zeros(N);
    w(1:Nfacet(1), 1:Nfacet(2)) = ones(Nfacet);

    % Generate the 2D Gaussian window
    alpha_x = 2*(N(2) - 1)*sigx;
    alpha_y = 2*(N(1) - 1)*sigy;
    wy = gausswin(N(1), alpha_y);
    wx = gausswin(N(2), alpha_x);
    wg = wy.*(wx.');

    % Convolve with a unit window of the facet size
    w = conv2(wg, w, 'same');

    % Threshold the smallest coefficients to zero / renormalize the window
    w(w < 1e-5) = 0; % set threshold as a parameter
    w = w/max(w(:)); % normalized window
end

end
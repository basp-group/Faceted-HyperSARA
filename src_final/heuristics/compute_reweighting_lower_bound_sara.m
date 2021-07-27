function [sig, mu] = ...
    compute_reweighting_lower_bound_sara(sigma_noise, squared_operator_norm)

% compute sig (estimate of the "noise level" in SARA space) involved in the
% reweighting scheme
sig = sqrt((sigma_noise^2)/squared_operator_norm);
mu = sig;

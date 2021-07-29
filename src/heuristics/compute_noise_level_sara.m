function sig = compute_noise_level_sara(sigma_noise, ...
    squared_operator_norm)
% Estimate noise level in the SARA domain.
%
% Return an estimate of the noise level transferred successively from the
% data to the image domain, then to the SARA domain.
%
% Parameters
% ----------
% sigma_noise : double
%     Noise variance in the data domain.
% squared_operator_norm : double
%     Squared norm of the measurement operator.
%
% Returns
% -------
% sig : double
%     Estimate of the noise level.
%

sig = sqrt((sigma_noise^2)/squared_operator_norm);

end

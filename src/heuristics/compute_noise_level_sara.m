function sig = compute_noise_level_sara(std_noise, squared_operator_norm)
% Estimate noise level in the SARA domain.
%
% Return an estimate of the noise level (std) transferred successively from
% the data to the image domain, then to the SARA domain.
%
% Parameters
% ----------
% std_noise : double
%     Noise standrad deviation in the data domain.
% squared_operator_norm : double
%     Squared norm of the measurement operator.
%
% Returns
% -------
% sig : double
%     Estimate of the noise level.
%

sig = sqrt((std_noise^2)/squared_operator_norm);

end

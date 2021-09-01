function [v0, weights0] = initialize_nuclear_serial(x, reweighting_alpha, sig_bar)
% Initalize the dual variables related to the l21-norm prior.
% 
% Parameters
% ----------
% x : array (3d)
%     Wideband image [M, N, L].
% reweighting_alpha : double
%     Reweighting parameter.
% sig_bar : double
%     Reweighting floor level.
%
% Returns
% -------
% v0 : cell
%     Dual variable associated with the nuclear-norm  prior 
%     [min(M*N, L), 1].
% weights0 : cell 
%     Associated weights for the reweigthing step.
%
              
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [26/03/2021]
%-------------------------------------------------------------------------%
%%

[M, N, c] = size(x);

v0 = zeros(M*N, c);

x0 = reshape(x, [M*N, c]);
[~,S0,~] = svd(x0,'econ');
d0 = abs(diag(S0));
upsilon_bar = sig_bar*reweighting_alpha;
weights0 = upsilon_bar ./ (upsilon_bar + d0);
    
end
        
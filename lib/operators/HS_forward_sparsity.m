function x = HS_forward_sparsity(x_wave, Psi, n, m)
% [extended_summary]
% 
% Parameters
% ----------
% x_wave : [type]
%     [description]
% Psi : [type]
%     [description]
% n : [type]
%     [description]
% m : [type]
%     [description]
% 
% Returns
% -------
% [type]
%     [description]

    [~, c] = size(x_wave);
    x = zeros(n, m, c);

    for i = 1:c
        %     dwtmode('per');
        x(:, :, i) = Psi(x_wave(:, i));
    end

end

function x_wave = HS_adjoint_sparsity(x, Psit, b)
% [extended_summary]
% 
% Parameters
% ----------
% x : [type]
%     [description]
% Psit : [type]
%     [description]
% b : [type]
%     [description]
% 
% Returns
% -------
% [type]
%     [description]

    [~, ~, c] = size(x);
    x_wave = zeros(b, c);

    for i = 1:c
        x_wave(:, i) = Psit(x(:, :, i));
    end

end

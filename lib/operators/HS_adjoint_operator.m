function x = HS_adjoint_operator(y, Gw, At, N, M)
% [extended_summary]
% 
% Parameters
% ----------
% y : [type]
%     [description]
% Gw : [type]
%     [description]
% At : [type]
%     [description]
% N : [type]
%     [description]
% M : [type]
%     [description]
% 
% Returns
% -------
% [type]
%     [description]

    c = length(y);
    x = zeros(N, M, c);

    for ind = 1:c
        x(:, :, ind) = At(Gw{ind}' * y{ind});
    end

end

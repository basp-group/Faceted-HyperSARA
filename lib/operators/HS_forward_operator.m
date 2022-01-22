function y = HS_forward_operator(x, Gw, A)
% [extended_summary]
% 
% Parameters
% ----------
% x : [type]
%     [description]
% Gw : [type]
%     [description]
% A : [type]
%     [description]
% 
% Returns
% -------
% [type]
%     [description]

    [~, ~, c] = size(x);
    y = cell(c, 1);

    for ind = 1:c
        y{ind} =  Gw{ind} * A(x(:, :, ind));
    end

end

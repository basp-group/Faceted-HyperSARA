function y = HS_forward_operator_precond(x, Gw, A, aW)
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
% aW : [type]
%     [description]
% 
% Returns
% -------
% [type]
%     [description]

    [~, ~, c] = size(x);
    y = cell(c, 1);

    for ind = 1:c
        y{ind} =  sqrt(cell2mat(aW{ind})) .* (Gw{ind} * A(x(:, :, ind)));
    end

end

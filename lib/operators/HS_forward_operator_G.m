function y = HS_forward_operator_G(x, G, W, A,flag_dr,Sigma)
% [extended_summary]
% 
% Parameters
% ----------
% x : [type]
%     [description]
% G : [type]
%     [description]
% W : [type]
%     [description]
% A : [type]
%     [description]
% flag_dr : [type]
%     [description]
% Sigma : [type]
%     [description]
% 
% Returns
% -------
% [type]
%     [description]

if nargin ==4
    flag_dr=0;
    Sigma =[];
elseif nargin ==5,Sigma =[];
end
[~, ~, c] = size(x);
y = cell(c, 1);

for i = 1:c
    Fx = A(x(:, :, i));
    for j = 1:length(G{i})
        if flag_dr
            if istril(G{i}{j})
                FxSlice = Fx(W{i}{j});
                y{i}{j} =Sigma{i}{j}.*( G{i}{j} * FxSlice + ( FxSlice' * G{i}{j})') ; FxSlice =[];
            else  y{i}{j} = Sigma{i}{j}.*( G{i}{j} * Fx(W{i}{j}));
            end
        else,  y{i}{j} = (G{i}{j} * Fx(W{i}{j}));
        end
    end
end

end

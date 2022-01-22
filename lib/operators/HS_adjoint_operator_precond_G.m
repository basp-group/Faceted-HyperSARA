function x = HS_adjoint_operator_precond_G(y, G, W, At, aW, N, M,flag_dr,Sigma)
% [extended_summary]
% 
% Parameters
% ----------
% y : [type]
%     [description]
% G : [type]
%     [description]
% W : [type]
%     [description]
% At : [type]
%     [description]
% aW : [type]
%     [description]
% N : [type]
%     [description]
% M : [type]
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

if nargin ==7
    flag_dr=0;
    Sigma =[];
elseif nargin ==8,Sigma =[];
end
c = length(y);
x = zeros(N, M, c);
% No = size(G{1}{1}, 2);
No = size(W{1}{1}, 1);

for i = 1:c
    g2 = zeros(No, 1);
    for j = 1:length(G{i})
        if flag_dr
            if istril(G{i}{j})
                weighted_y = (sqrt(aW{i}{j}).*Sigma{i}{j}) .* y{i}{j};
                g2(W{i}{j}) = g2(W{i}{j}) + G{i}{j}' * weighted_y  +  G{i}{j} * weighted_y;
                
            else, g2(W{i}{j}) = g2(W{i}{j}) + G{i}{j}' * ((sqrt(aW{i}{j}).*Sigma{i}{j}) .* y{i}{j}) ;%
            end
        else
            g2(W{i}{j}) = g2(W{i}{j}) + G{i}{j}' * (sqrt(aW{i}{j}) .* y{i}{j});
        end
    end
    x(:, :, i) = At(g2);
end

end

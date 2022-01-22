function x = HS_adjoint_operator_G(y, G, W, At, N, M,flag_dr,Sigma)
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


if nargin ==6
    flag_dr=0;
    Sigma =[];
elseif nargin ==7,Sigma =[];
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
		     weighted_data = (Sigma{i}{j}.*y{i}{j});
		     g2(W{i}{j}) = g2(W{i}{j}) + G{i}{j}' * weighted_data  +  G{i}{j} * weighted_data;
	       else, g2(W{i}{j}) = g2(W{i}{j}) + G{i}{j}' * (Sigma{i}{j}.*y{i}{j});
               end
	    else, g2(W{i}{j}) = g2(W{i}{j}) + G{i}{j}' * (y{i}{j});
            end
        end
        x(:, :, i) = At(g2);
    end

end

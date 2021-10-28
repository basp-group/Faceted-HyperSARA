function y = HS_forward_operator_G(x, G, W, A,flag_dr)
if nargin <5
    flag_dr=0;
end
    [~, ~, c] = size(x);
    y = cell(c, 1);

    for i = 1:c
        Fx = A(x(:, :, i));
        for j = 1:length(G{i})
            %         y{i}{j} = sqrt(aW{i}{j}) .* (G{i}{j} * Fx);
            if flag_dr
                 y{i}{j} = G{i}{j} * Fx(W{i}{j}) + G{i}{j}' * Fx(W{i}{j});
            else,  y{i}{j} = (G{i}{j} * Fx(W{i}{j}));
            end
        end
    end

end

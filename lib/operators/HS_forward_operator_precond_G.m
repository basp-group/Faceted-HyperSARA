function y = HS_forward_operator_precond_G(x, G, W, A, aW,flag_dr)
if nargin <6
    flag_dr=0;
end
    [~, ~, c] = size(x);
    y = cell(c, 1);

    for i = 1:c
        Fx = A(x(:, :, i));

        y{i}  = cell(size(G{i}));
        for j = 1:length(G{i})
            if flag_dr
                y{i}{j} = sqrt(aW{i}{j}) .* (G{i}{j} * Fx(W{i}{j}) + G{i}{j}' * Fx(W{i}{j}));
            else
            %         y{i}{j} = sqrt(aW{i}{j}) .* (G{i}{j} * Fx);
            y{i}{j} = sqrt(aW{i}{j}) .* (G{i}{j} * Fx(W{i}{j}));
            end
        end
    end

end

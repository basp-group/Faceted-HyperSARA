function x = HS_adjoint_operator_precond_G(y, G, W, At, aW, N, M,flag_dr)
if nargin <8
    flag_dr=0;
end
    c = length(y);
    x = zeros(N, M, c);
    % No = size(G{1}{1}, 2);
    No = size(W{1}{1}, 1);

    for i = 1:c
        g2 = zeros(No, 1);
        for j = 1:length(G{i})
            if flag_dr
                g2(W{i}{j}) = g2(W{i}{j}) + G{i}{j}' * (sqrt(aW{i}{j}) .* y{i}{j}) + G{i}{j} * (sqrt(aW{i}{j}) .* y{i}{j});
            else
            %         g2 = g2 + G{i}{j}' * (sqrt(aW{i}{j}) .* y{i}{j});
                g2(W{i}{j}) = g2(W{i}{j}) + G{i}{j}' * (sqrt(aW{i}{j}) .* y{i}{j});
            end
        end
        x(:, :, i) = At(g2);
    end

end

function y = HS_forward_operator(x,Gw,A)

% Parameters
[~, ~, c] = size(x);
y = zeros(size(Gw{1},1), c);

%
for ind = 1:c
    y(:, ind) = Gw{ind} * A(x(:,:,ind));
end

end

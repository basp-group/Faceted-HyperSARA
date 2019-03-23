function y = HS_fouRed_adjoint_operator(x, Bt, aW)

% Parameters
[~, c] = size(x);

%
for ind = 1:c
    if exist('aW', 'var')
        y(:,:,ind) = Bt{ind}(sqrt(cell2mat(aW{ind})) .* x(:,ind));
    else
        y(:,:,ind) = Bt{ind}(x(:,ind));
    end
end

end

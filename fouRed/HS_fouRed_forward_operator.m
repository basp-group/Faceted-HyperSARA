function y = HS_fouRed_forward_operator(x, B, aW)

% Parameters
[~, ~, c] = size(x);

%
for ind = 1:c
    if exist('aW', 'var')
        y(:,ind) =  sqrt(cell2mat(aW{ind})) .* B{ind}(x(:,:,ind));
    else
        y(:,ind) =  B{ind}(x(:,:,ind));
    end
end

end

function x = HS_forward_sparsity(y,Psi,n,m)

% Parameters
[~,c] = size(y);
x = zeros(n,m,c);

%
for i = 1 : c
    %     dwtmode('per');
    x(:,:,i) = Psi(y(:,i));
end

end

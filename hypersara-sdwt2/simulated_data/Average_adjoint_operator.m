function x = Average_adjoint_operator(x_ave,cent,step,ave,n,m,c)

%
x = zeros(n,m,c);

for i = 1 : length(cent)
    x(:,:,cent(i)-step:cent(i)+step) = repmat(x_ave(:,:,i),1,1,ave);
end

end
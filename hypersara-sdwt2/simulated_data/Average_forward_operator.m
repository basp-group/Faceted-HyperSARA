function x_ave = Average_forward_operator(x,cent,step,n,m)

%
x_ave = zeros(n,m,length(cent));

for i = 1 : length(cent)
    x_ave(:,:,i) = mean(x(:,:,cent(i)-step:cent(i)+step),3);
end

end
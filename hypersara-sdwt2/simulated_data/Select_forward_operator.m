function x_ave = Select_forward_operator(x,cent,n,m)

%
x_ave = zeros(n,m,length(cent));

for i = 1 : length(cent)
    x_ave(:,:,i) = x(:,:,cent(i));
end

end
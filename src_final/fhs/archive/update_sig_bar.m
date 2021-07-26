function sig_bar = update_sig_bar(x, b)

nChannels = size(x, 3);
b = reshape(b, [size(x,1)*size(x,2), nChannels]);
[U,~,V] = svd(reshape(x, [size(x,1)*size(x,2), nChannels]),'econ');
sig_bar = std(abs(diag(U'*b*V)));

end
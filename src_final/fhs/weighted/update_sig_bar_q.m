function sig_bar_q = update_sig_bar_q(xq, wq, bq)

nChannels = size(xq, 3);
bq = reshape(bq, [size(xq,1)*size(xq,2), nChannels]);
[U,~,V] = svd(reshape(xq.*wq, [size(xq,1)*size(xq,2), nChannels]),'econ');
sig_bar_q = std(abs(diag(U'*bq*V)));

end
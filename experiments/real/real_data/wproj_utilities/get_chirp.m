function wkernel =get_chirp(W,nTerm,PHASE,levelC)
chirp  = exp(-2*1i*pi*W.*nTerm); % chirp analytical expr. in the image space
chirpFourier = ifftshift(PHASE.*fft2(chirp)); % no phase is adopted in the FFT (check NUFFT shift convention)
chirpFourier =  chirpFourier./norm(chirpFourier(:));
% sparsification of the chirp
[thresVal,~] =  BISECT(levelC, 1e-5, abs(chirpFourier),2);
wkernel = full(chirpFourier.*(abs(chirpFourier)>thresVal));
wkernel = wkernel./(norm(nonzeros(wkernel))); % normalisation to avoid energy loss (not compulsory)
end
function normalise_residual(filename, outputname, psfDir, ch)
fprintf("channels: %d\n", ch)
nChannels = length(ch);
res = fitsread(filename);
res_norm = zeros(size(res));
fprintf("Read file: %s\n\n", filename)
for i = 1:nChannels
    psfName = [psfDir, '/psf_2b_ind1_16=', num2str(ch(i)), '.fits'];
    fprintf("Read psf: %s\n\n", psfName)
    psf = fitsread(psfName);
    peak = max(abs(psf(:)));
    res_norm(:,:,i) = res(:,:,i) / peak;
end
fitswrite(res_norm, outputname)
fprintf("Normalisation of psf is finished!\n\n")
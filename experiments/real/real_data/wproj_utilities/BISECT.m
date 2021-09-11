function [threshold, n] = BISECT(E_threshold, tolerance, abschirp, L_norm)
% calculates the threshold for keeping E_hreshold= percentile  of total
% energy of the chirp kernel, where abschirp is the absolute value of the
% row of the convoution kernel
% L_norm =1;
x0 = 0.0;
x1 = 1;
maxval = max(max(abschirp));
n = 1;
if L_norm == 2; Energy0 = sqrt(sum(sum(abschirp.^2)));
else;  Energy0 = (sum(sum(abschirp)));
end

Energy = 1;
while (abs(x0 - x1) / abs(x1)) >= tolerance || (Energy - E_threshold) < 0.0
    x = x0 + (x1 - x0) / 2.0;
    thresvalue = x * maxval;
    if L_norm == 2; Energy = sqrt(sum(sum(abschirp(abschirp > thresvalue).^2))) ./ Energy0;
    else;   Energy = (sum(sum(abschirp(abschirp > thresvalue)))) ./ Energy0;
    end
    if Energy < E_threshold;      x1 = x;
    else; x0 = x;
    end
    xtest = x0 + (x1 - x0) / 2;
    thresvalue = xtest * maxval;
    if L_norm == 2;   Energy = sqrt(sum(sum(abschirp(abschirp > thresvalue).^2))) ./ Energy0;
    else;  Energy = (sum(sum(abschirp(abschirp > thresvalue)))) ./ Energy0;
    end
    n = n + 1;
    if n > 1000;   break
    end
end
threshold = thresvalue;
end

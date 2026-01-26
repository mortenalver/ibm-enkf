
% Test how the periodograms look like for series made up of periodcally
% independent values

nparts = 1000;
nper = 15;

t1 = zeros(nparts*nper,1);

for i=1:nparts

    meanval = rand(1);
    t1((1+(i-1)*nper):i*nper) = meanval + 0.015*randn(nper,1);

end

%figure, plot(t1)

[psd_t1, w1] = periodogram(t1, hann(length(t1)));

figure,semilogy(w1, psd_t1)
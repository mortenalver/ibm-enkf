
ntime = 1000;
N = 5;
sigma = 1;
a = zeros(ntime,N);
a(1,:) = sigma*randn(1,N);
for i=2:ntime
    a(i,:) = a(i-1,:) + sigma*randn(1,N);
end

hw = hann(ntime);
period =[];
for i=1:N
    [pxx, w] = periodogram(a(:,i),hw);
    if i==1
        period = pxx;
        %figure,semilogy(w, pxx);
    else
        period = period + pxx;
    end
end
period = period/N;

figure,semilogy(w, period)

clc, clear;
N=1e6;
SNR=-4:2:14;
snr=db2pow(SNR);
BER=zeros([1 length(SNR)]);
ber=zeros([1 length(SNR)]);

for u=1:length(snr)
    a = randi([0 1],1,N);
    x = 2*a-1;
    %AWGN standard
    w=(1/sqrt(2*snr(u)))*randn(1,N);
    y1 = x+w;
    %AWGN function
    y2=awgn(x,SNR(u),'measured');
    
    x1_hat = ones(1,N);         
    x1_hat(y1 < 0) = 0;
    
    x2_hat = ones(1,N);         
    x2_hat(y2 < 0) = 0;
    
    BER(u) = biterr(a,x1_hat)/N;
    ber(u) = biterr(a,x2_hat)/N;
end
figure (1)
colorstring = 'rk';%ymcrgbwk
semilogy(SNR, BER,'s-','Color', colorstring(1), 'linewidth', 1.5);
hold on; grid on;
semilogy(SNR, ber,'d--','Color', colorstring(2), 'linewidth', 1.5); 
legend('Standard','AWGN function');
title('Standard AWGN vs AWGN function(BPSK BER)');
xlabel('EbNo(dB)')   
ylabel('BER')
grid on
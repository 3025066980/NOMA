clc; clear variables;
N = 10^5;         
SNR = 0:40;      
P=1;
a1 = 0.75; a2 = 0.18; a3 = 0.05; a4 = 0.02;
ber1=zeros([1 length(SNR)]);
ber2=zeros([1 length(SNR)]);
ber3=zeros([1 length(SNR)]);
ber4=zeros([1 length(SNR)]);
for u = 1:length(SNR)
    x1 = randi([0 1],1,N);
    x2 = randi([0 1],1,N);
    x3 = randi([0 1],1,N);
    x4 = randi([0 1],1,N);    
    
    xmod1 = 2*x1-1;
    xmod2 = 2*x2-1;
    xmod3 = 2*x3-1;
    xmod4 = 2*x4-1;
    x = sqrt(a1*P)*xmod1 + sqrt(a2*P)*xmod2+ sqrt(a3*P)*xmod3 + sqrt(a4*P)*xmod4;

    y1 = awgn(x,SNR(u),'measured');
    y2 = awgn(x,SNR(u),'measured');
    y3 = awgn(x,SNR(u),'measured');
    y4 = awgn(x,SNR(u),'measured');

    %AT USER 1
    %Direct decoding of x from y1
    x1_hat = ones(1,N);         
    x1_hat(y1 < 0) = 0;         %Final bits for user 1
    
    %AT USER 2
    rem2=y2-sqrt(a1*P)*xmod1;
    
    %Decode x2 from rem without signal x1
    x2_hat = zeros(1,N);
    x2_hat(rem2>0) = 1;          %Final bits for user 2
    
    %AT USER 3
    rem3 = y3 - sqrt(a1*P)*xmod1 - sqrt(a2*P)*xmod2;
    
    %Decode x3 from rem
    x3_hat = zeros(1,N);
    x3_hat(rem3>0) = 1;          %Final bits for user 3
    
    %AT USER 4
    %Direct decoding of x from y2(dd3) First of all, estimate the previous
    rem4=y4 - sqrt(a1*P)*xmod1 - sqrt(a2*P)*xmod2 - sqrt(a3*P)*xmod3;
    
    x4_hat = zeros(1,N);
    x4_hat(rem4>0) = 1;
    
    %Estimate BER
    ber1(u) = biterr(x1,x1_hat)/N;
    ber2(u) = biterr(x2,x2_hat)/N;
    ber3(u) = biterr(x3,x3_hat)/N;
    ber4(u) = biterr(x4,x4_hat)/N;
end
%plot BER curves
colorstring = 'bmrg';
figure(1)
semilogy(SNR, ber1,'+--','Color', colorstring(1), 'linewidth', 1); 
hold on; grid on;
semilogy(SNR, ber2,'+--','Color', colorstring(2), 'linewidth', 1);
semilogy(SNR, ber3,'+--','Color', colorstring(3), 'linewidth', 1); 
semilogy(SNR, ber4,'+--','Color', colorstring(4), 'linewidth', 1); 
title('BER graph for NOMA in AWGN channel-perfect SIC-experimental');
xlabel('SNR');
ylabel('BER');
clc; clear variables;
N = 10^5; SNR = 0:40; P=1;
a1 = 0.70; a2 = 0.25; a3 = 0.05;
ber1=zeros([1 length(SNR)]);
ber2=zeros([1 length(SNR)]);
ber3=zeros([1 length(SNR)]);
a=zeros(1,length(SNR));
% b=zeros(1,length(SNR));
c=zeros(1,length(SNR));
for u = 1:length(SNR)
    x1 = randi([0 1],1,N);
    x2 = randi([0 1],1,N);
    x3 = randi([0 1],1,N);
    xmod1 = 2*x1-1;
    xmod2 = 2*x2-1;
    xmod3 = 2*x3-1;
    xx = sqrt(a1*P)*xmod1 + sqrt(a2*P)*xmod2+ sqrt(a3*P)*xmod3;
    xx3=sqrt(a3*P)*xmod3;
    noise1 = wgn(1,length(xx),10*log10(0.7/(10^(SNR(u)/10))));
    noise2 = wgn(1,length(xx),10*log10(0.25/(10^(SNR(u)/10))));
    noise3 = wgn(1,length(xx),10*log10(0.05/(10^(SNR(u)/10))));%0.5

    y11 = xx + noise1;
    y22 = xx + noise2;
    y33 = xx + noise3;
    y331=xx3 + noise3;

    %AT USER 1
    %Direct decoding of x from y1
    x1_hat = ones(1,N);         
    x1_hat(y11 < 0) = 0;%Final bits for user 1
    
    %AT USER 2
    x21_est = ones(1,N);
    x21_est(y22 < 0) = -1;       %Estimate user 1's signal first
    rem2 = y22 - sqrt(a1*P)*x21_est;%Extract useless signal 1 from y2
    %Decode x2 from rem without signal x1
    x2_hat = zeros(1,N);
    x2_hat(rem2>0) = 1;%Final bits for user 2
    
    %AT USER 3
    x31_est = ones(1,N);
    x31_est(y33 < 0) = -1; %Estimate user 1's signal first
    rem23 = y33 - sqrt(a1*P)*x31_est; %Extract useless signal 1 from y3
    a(u)=mean(xmod1-x31_est*sqrt(a1*P));
    x32_est = ones(1,N);
    x32_est(rem23<0) = -1;
    rem33=rem23 - sqrt(a2*P)*x32_est; %Extract useless signal 2 from y3
    b(u,:)=xmod3-rem33;
    %Decode x3 from rem
    x3_hat = zeros(1,N);
    x3_hat(rem33>0) = 1;%Final bits for user 3

    x43_hat = ones(1,N);         
    x43_hat(y331 < 0) = 0;%Final bits for user 1
    %Estimate BER
    ber1(u) = biterr(x1,x1_hat)/N;
    ber2(u) = biterr(x2,x2_hat)/N;
    ber3(u) = biterr(x3,x3_hat)/N;
    ber31(u) = biterr(x3,x43_hat)/N;
end
% plot BER curves
colorstring = 'bmry';
figure(2)
semilogy(SNR, ber1,'+--','Color', colorstring(1), 'linewidth', 1); 
hold on; grid on;
semilogy(SNR, ber2,'+--','Color', colorstring(2), 'linewidth', 1);
semilogy(SNR, ber3,'+--','Color', colorstring(3), 'linewidth', 1); 
semilogy(SNR, ber31,'+--','Color', colorstring(4), 'linewidth', 1); 
legend('User 1','User 2','User 3','Location','northeast');
title('BER graph for NOMA in AWGN channel-imperfect SIC-experimental');
xlabel('SNR(dB)');
ylabel('BER');
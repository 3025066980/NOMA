clc;clear;
N = 10^5;         
SNR = 0:50;      
P=1;
a1 = 0.70; a2 = 0.25; a3 = 0.05;
j=100;
bber1=zeros(j,51);
bber2=zeros(j,51);
bber3=zeros(j,51);
for kk = 1:j
    h = raylrnd(0.75,1,3);
    h = sort(h);

    for u = 1:length(SNR)
        x1 = randi([0 1],1,N);
        x2 = randi([0 1],1,N);
        x3 = randi([0 1],1,N); 
    
        xmod1 = 2*x1-1;
        xmod2 = 2*x2-1;
        xmod3 = 2*x3-1;

        xx1 = h(1)*(sqrt(a1*P)*xmod1 + sqrt(a2*P)*xmod2+ sqrt(a3*P)*xmod3);
        xx2 = h(2)*(sqrt(a1*P)*xmod1 + sqrt(a2*P)*xmod2+ sqrt(a3*P)*xmod3);
        xx3 = h(3)*(sqrt(a1*P)*xmod1 + sqrt(a2*P)*xmod2+ sqrt(a3*P)*xmod3);
        
        noise = wgn(1,length(xx1),10*log10(1/(10^(SNR(u)/10))));

        y11 = xx1 + noise;
        y22 = xx2 + noise;
        y33 = xx3 + noise;    
        y1 = y11/h(1);
        y2 = y22/h(2);
        y3 = y33/h(3);

        %AT USER 1
        %Direct decoding of x from y1
        x1_hat = ones(1,N);
        x1_hat(y1 < 0) = 0;%Final bits for user 1

        %AT USER 2
        x21_est = ones(1,N);
        x21_est(y2 < 0) = -1;       %Estimate user 1's signal first
        rem2 = y2 - sqrt(a1*P)*x21_est;%Extract useless signal 1 from y2
        %Decode x2 from rem without signal x1
        x2_hat = zeros(1,N);
        x2_hat(rem2>0) = 1;%Final bits for user 2

        %AT USER 3
        x31_est = ones(1,N);
        x31_est(y3 < 0) = -1; %Estimate user 1's signal first
        rem23 = y3 - sqrt(a1*P)*x31_est; %Extract useless signal 1 from y3
        x32_est = ones(1,N);
        x32_est(rem23<0) = -1;
        rem33=rem23 - sqrt(a2*P)*x32_est; %Extract useless signal 2 from y3
        %Decode x3 from rem
        x3_hat = zeros(1,N);
        x3_hat(rem33>0) = 1;%Final bits for user 3

        %Estimate BER
        bber1(kk,u) = biterr(x1,x1_hat)/N;
        bber2(kk,u) = biterr(x2,x2_hat)/N;
        bber3(kk,u) = biterr(x3,x3_hat)/N;
    end
end
ber1 = mean(bber1);
ber2 = mean(bber2);
ber3 = mean(bber3);
colorstring = 'bmr';
figure(4)
semilogy(SNR, ber1,'+--','Color', colorstring(1), 'linewidth', 1); 
hold on; grid on; ylim([10^(-5) 1]);
semilogy(SNR, ber2,'+--','Color', colorstring(2), 'linewidth', 1);
semilogy(SNR, ber3,'+--','Color', colorstring(3), 'linewidth', 1); 
title('BER graph for NOMA in AWGN channel-perfect SIC-experimental');
xlabel('SNR');
ylabel('BER');
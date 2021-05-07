clc;
clear;
N = 10^5;         
SNR = 0:40;      
P=1;
a1 = 0.75; a2 = 0.18; a3 = 0.05; a4 = 0.02;
bber1=zeros(100,41);
bber2=zeros(100,41);
bber3=zeros(100,41);
bber4=zeros(100,41);
for kk = 1:100
    h = raylrnd(0.75,1,4);
    h = sort(h);

    for u = 1:length(SNR)
        x1 = randi([0 1],1,N);
        x2 = randi([0 1],1,N);
        x3 = randi([0 1],1,N);
        x4 = randi([0 1],1,N);    
    
        xmod1 = 2*x1-1;
        xmod2 = 2*x2-1;
        xmod3 = 2*x3-1;
        xmod4 = 2*x4-1;

        xx1 = h(1)*(sqrt(a1*P)*xmod1 + sqrt(a2*P)*xmod2+ sqrt(a3*P)*xmod3 + sqrt(a4*P)*xmod4);
        xx2 = h(2)*(sqrt(a1*P)*xmod1 + sqrt(a2*P)*xmod2+ sqrt(a3*P)*xmod3 + sqrt(a4*P)*xmod4);
        xx3 = h(3)*(sqrt(a1*P)*xmod1 + sqrt(a2*P)*xmod2+ sqrt(a3*P)*xmod3 + sqrt(a4*P)*xmod4);
        xx4 = h(4)*(sqrt(a1*P)*xmod1 + sqrt(a2*P)*xmod2+ sqrt(a3*P)*xmod3 + sqrt(a4*P)*xmod4);
        
        noise = wgn(1,length(xx1),10*log10(1/(10^(SNR(u)/10))));

        y11 = xx1 + noise;
        y22 = xx2 + noise;
        y33 = xx3 + noise;
        y44 = xx4 + noise;
    
        y1 = y11/h(1);
        y2 = y22/h(2);
        y3 = y33/h(3);
        y4 = y44/h(4);

        %AT USER 1
        %Direct decoding of x from y1
        x1_hat = y1>0;        

        %AT USER 2
        rem2=y2 - sqrt(a1*P)*xmod1;

        %Decode x2 from rem without signal x1
        x2_hat = rem2>0;

        %AT USER 3
        rem3 = y3 - sqrt(a1*P)*xmod1 - sqrt(a2*P)*xmod2;

        %Decode x3 from rem
        x3_hat = rem3>0;

        %AT USER 4
        %Direct decoding of x from y2(dd3) First of all, estimate the previous
        rem4=y4 - sqrt(a1*P)*xmod1 - sqrt(a2*P)*xmod2 - sqrt(a3*P)*xmod3;
        x4_hat = rem4>0;

        %Estimate BER
        bber1(kk,u) = biterr(x1,x1_hat)/N;
        bber2(kk,u) = biterr(x2,x2_hat)/N;
        bber3(kk,u) = biterr(x3,x3_hat)/N;
        bber4(kk,u) = biterr(x4,x4_hat)/N;
    end
end

ber1 = mean(bber1);
ber2 = mean(bber2);
ber3 = mean(bber3);
ber4 = mean(bber4);

colorstring = 'bmrg';
figure(2)
semilogy(SNR, ber1,'+--','Color', colorstring(1), 'linewidth', 1); 
hold on; grid on;
semilogy(SNR, ber2,'+--','Color', colorstring(2), 'linewidth', 1);
semilogy(SNR, ber3,'+--','Color', colorstring(3), 'linewidth', 1); 
semilogy(SNR, ber4,'+--','Color', colorstring(4), 'linewidth', 1); 
title('BER graph for NOMA in AWGN channel-perfect SIC-experimental');
xlabel('SNR');
ylabel('BER');
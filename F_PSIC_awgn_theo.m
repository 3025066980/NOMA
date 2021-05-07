clear, clc;
P=1;
[a1, a2, a3, a4] = deal(0.75 , 0.18 , 0.05 , 0.02);
[A1, A2 ,A3 ,A4] = deal(sqrt(P*a1),sqrt(P*a2),sqrt(P*a3),sqrt(P*a4));
A=[A1, A2 ,A3 ,A4];
SNR = 0:40;             %SNR range in dB
variance=10.^(-SNR./10);
mtx1=[1 1 1 1;1 -1 1 1;1 -1 -1 1;1 1 -1 1; 1 -1 1 -1;1 1 1 -1; 1 1 -1 -1;1 -1 -1 -1];
mtx2=[1 1 1;1 -1 1;1 1 -1;1 -1 -1];
mtx3=[1 1;1 -1];
mtx4=1;
M=4;
[ls1, ls2, ls3, ls4]=deal(8, 4, 2, 1);
[sumt1, sumt2,sumt3,sumt4]=deal(0,0,0,0);
gamma=P./variance;
%correct power
for h =1:M 
        mtx1(:,h)=mtx1(:,h)*A(h);
end
for h =1:M-1 
        mtx2(:,h)=mtx2(:,h)*A(h+1);
end
for h =1:M-2 
        mtx3(:,h)=mtx3(:,h)*A(h+2);
end
mtx4=A(4)*mtx4;
%Summation
for j=1:ls1
    Y=mtx1(j,:);
    beta=sum(Y)/sqrt(P);
    sumt1=sumt1+qfunc(beta*sqrt(gamma));
end
for j=1:ls2
    Y=mtx2(j,:);
    beta=sum(Y)/sqrt(P);
    sumt2=sumt2+qfunc(beta*sqrt(gamma));
end
for j=1:ls3
    Y=mtx3(j,:);
    beta=sum(Y)/sqrt(P);
    sumt3=sumt3+qfunc(beta*sqrt(gamma));
end
for j=1:ls4
    Y=mtx4(j,:);
    beta=sum(Y)/sqrt(P);
    sumt4=sumt4+qfunc(beta*sqrt(gamma));
end
k1=1/(ls1)*sumt1;
k2=1/(ls2)*sumt2;
k3=1/(ls3)*sumt3;
k4=1/(ls4)*sumt4;
colorstring = 'bmrgk';
figure(1)
semilogy(SNR,k1,'.-','Color', colorstring(1),'LineWidth',1)
ylim([10^(-5) 1])
hold on;grid on;
semilogy(SNR,k2,'.-','Color', colorstring(2),'LineWidth',1)
semilogy(SNR,k3,'.-','Color', colorstring(3),'LineWidth',1)
semilogy(SNR,k4,'.-','Color', colorstring(4),'LineWidth',1)
snr = db2pow(SNR);      %SNR range in linear scale
theoryBer = qfunc(sqrt(1*snr));
semilogy(SNR+17,theoryBer,'*-','Color',colorstring(5),'linewidth',1);
legend('User 1 \alpha_1 = 0.75','User 2 \alpha_2 = 0.18','User 3 \alpha_3 = 0.05','User 4 \alpha_4 = 0.02','BER Theoritical','Location','southwest');
txt1 = 'Simulation -----,Theoritical: Solid Line';
text(6,3.43e-05,txt1)
title('BER graph for NOMA in AWGN channel-theoretical');
xlabel('SNR');
ylabel('BER');
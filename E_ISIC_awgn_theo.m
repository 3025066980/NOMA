clear, clc;
P=1;
[a1, a2, a3] = deal(0.70, 0.25, 0.05);
a=[a1, a2, a3];
[A1, A2 ,A3] = deal(sqrt(P*a1),sqrt(P*a2),sqrt(P*a3));
A=[A1, A2 ,A3];
SNR = 0:50;
variance=10.^(-SNR./10);
mtx1=[1 1 1;1 -1 1;1 1 -1;1 -1 -1];
mtx2=[1 1;1 -1];
mtx3=1;
M=3;OMEGA=1;zeta=0.001;
[ls1, ls2, ls3]=deal(4, 2, 1);
[sumt1, sumt2,sumt3]=deal(0,0,0);
gamma=P./variance;
%correct power
for h =1:M 
        mtx1(:,h)=mtx1(:,h)*A(h);
end
for h =1:M-1 
        mtx2(:,h)=mtx2(:,h)*A(h+1);
end
mtx3=A(3)*mtx3;
%Summation
for j=1:ls1
    
    Y=mtx1(j,:);
    beta=sum(Y)/sqrt(P);
    omega2=OMEGA*zeta*0;
    betaplus=beta.*sqrt(P./(P+gamma*omega2));
    sumt1=sumt1+qfunc(betaplus.*sqrt(gamma));
    
end
for j=1:ls2
    Y=mtx2(j,:);
    beta=sum(Y)/sqrt(P);
    omega2=OMEGA*zeta*(0.70);
    betaplus=beta.*sqrt(P./(P+gamma*omega2));%*0.5
    sumt2=sumt2+qfunc(betaplus.*sqrt(gamma));

end
for j=1:ls3
    Y=mtx3(j,:);
    beta=sum(Y)/sqrt(P);
    omega2=OMEGA*zeta*(0.70+0.25);
    betaplus=beta.*sqrt(P./(P+gamma*omega2));%*0.7
    sumt3=sumt3+qfunc(betaplus.*sqrt(gamma));

end
k1=1/(ls1)*sumt1;
k2=1/(ls2)*sumt2;
k3=1/(ls3)*sumt3;
figure(1)
colorstring = 'bmr';
semilogy(SNR,k1,'-','Color', colorstring(1),'LineWidth',1)
ylim([10^(-5) 1])
hold on;grid on;
semilogy(SNR,k2,'-','Color', colorstring(2),'LineWidth',1)%+5
semilogy(SNR,k3,'-','Color', colorstring(3),'LineWidth',1)
legend('User 1 \alpha_1 = 0.70','User 2 \alpha_2 = 0.25','User 3 \alpha_4 = 0.05','Location','southwest');
txt1 = 'Simulation -----,Theoritical: Solid Line';
text(6,3.43e-05,txt1)
title('BER graph for NOMA in AWGN channel-theoretical');
xlabel('SNR');
ylabel('BER');
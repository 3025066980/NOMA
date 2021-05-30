clear, clc;
P=1;OMEGA=1;zeta=0.001;
[a1, a2, a3] = deal(0.70 , 0.25 , 0.05);
[A1, A2 ,A3] = deal(sqrt(P*a1),sqrt(P*a2),sqrt(P*a3));
A=[A1, A2 ,A3];
mtx1=[1 1 1;1 -1 1;1 1 -1;1 -1 -1];
mtx2=[1 1;1 -1];
mtx3=1;
omega21=OMEGA*zeta*0;
omega22=OMEGA*zeta*A1^2;
omega23=OMEGA*zeta*(A1^2+A2^2);
M=3;
[ls1, ls2, ls3]=deal(4, 2, 1);
[sumt1, sumt2,sumt3]=deal(0,0,0);
SNR = 0:40;
gamma=10.^(SNR./10);
variance=P./gamma;
gamma1=P./(variance + omega21);
gamma2=P./(variance + omega22);
gamma3=P./(variance + omega23);
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
    sumt1=sumt1+qfunc(beta*sqrt(gamma1));
end
for j=1:ls2
    Y=mtx2(j,:);
    beta=sum(Y)/sqrt(P);
    sumt2=sumt2+qfunc(beta*sqrt(gamma2));
end
for j=1:ls3
    Y=mtx3(j,:);
    beta=sum(Y)/sqrt(P);
    sumt3=sumt3+qfunc(beta*sqrt(gamma3));
end
k1=1/(ls1)*sumt1;
k2=1/(ls2)*sumt2;
k3=1/(ls3)*sumt3;
colorstring = 'bmr';
figure(2)
semilogy(pow2db(gamma1),k1,':','Color', colorstring(1),'LineWidth',1)
ylim([10^(-5) 1]);xlim([0 inf])
hold on;grid on;
semilogy(pow2db(gamma2),k2,':','Color', colorstring(2),'LineWidth',1)
hold on
semilogy(pow2db(gamma3),k3,':','Color', colorstring(3),'LineWidth',1)
hold on
legend('User 1 \alpha_1 = 0.70','User 2 \alpha_2 = 0.20','User 3 \alpha_3 = 0.05','Location','southwest');
txt1 = 'Simulation -----,Theoritical: Solid Line';
text(6,3.43e-05,txt1)
title('BER graph for NOMA in AWGN channel-theoretical');
xlabel('SNR');
ylabel('BER');
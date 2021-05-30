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
M=3;OMEGA=1;zeta=0;
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
vector1=zeros([1 length(SNR)]);
vector2=zeros([1 length(SNR)]);
vector3=zeros([1 length(SNR)]);
for h=1:length(SNR)
    sum1_3=0;
    m=1;
    for i=1:ls1
        sum1_2=0;
        for t=0:(M-m)
           sum1_1=0;
           for k=0:(m+t-1)
               ita=factorial(M)/(factorial(M-m)*factorial(m-1));
               Y=mtx1(i,:);
               beta=sum(Y)/sqrt(P);
               omega2=OMEGA*zeta*0;
               betaplus=beta.*sqrt(P./(P+gamma*omega2));
               omega=sqrt( 1+( (2*(k+1)) ./ (betaplus.^2*gamma(h)*1) ) );
               p1=(ita*(omega-1)*(-1)^(t+k)*nchoosek((M-m),t)*nchoosek((m+t-1),k)) / (2^(M-m+1)*(k+1)*omega);
               sum1_1=sum1_1+p1;
           end
           sum1_2=sum1_2+sum1_1;
        end
        sum1_3=sum1_3+sum1_2;
    end   
    vector1(h)=sum1_3;
end
for h=1:length(SNR)
    sum1_3=0;
    m=2;
    for i=1:ls2
        sum1_2=0;
        for t=0:(M-m)
           sum1_1=0;
           for k=0:(m+t-1)
               ita=factorial(M)/(factorial(M-m)*factorial(m-1));
               Y=mtx1(i,:);
               beta=sum(Y)/sqrt(P);
               omega2=OMEGA*zeta*(0.70);
               betaplus=beta.*sqrt(P./(P+gamma*omega2));%*0.5
               omega=sqrt( 1+( (2*(k+1)) ./ (betaplus.^2*gamma(h)*1) ) );
               p1=(ita*(omega-1)*(-1)^(t+k)*nchoosek((M-m),t)*nchoosek((m+t-1),k)) / (2^(M-m+1)*(k+1)*omega);
               sum1_1=sum1_1+p1;
           end
           sum1_2=sum1_2+sum1_1;
        end
        sum1_3=sum1_3+sum1_2;
    end   
    vector2(h)=sum1_3;
end
for h=1:length(SNR)
    sum1_3=0;
    m=3;
    for i=1:ls3
        sum1_2=0;
        for t=0:(M-m)
           sum1_1=0;
           for k=0:(m+t-1)
               ita=factorial(M)/(factorial(M-m)*factorial(m-1));
               Y=mtx2(i,:);
               beta=sum(Y)/sqrt(P);
               omega2=OMEGA*zeta*(0.70+0.25);
               betaplus=beta.*sqrt(P./(P+gamma*omega2));%*0.7
               omega=sqrt( 1+( (2*(k+1)) ./ (betaplus.^2*gamma(h)*1) ) );
               p1=(ita*(omega-1)*(-1)^(t+k)*nchoosek((M-m),t)*nchoosek((m+t-1),k)) / (2^(M-m+1)*(k+1)*omega);
               sum1_1=sum1_1+p1;
           end
           sum1_2=sum1_2+sum1_1;
        end
        sum1_3=sum1_3+sum1_2;
    end   
    vector3(h)=sum1_3;
end
colorstring = 'bmr';
figure(4)
b=(SNR+7.8)./SNR;
c=(SNR+11)./SNR;
semilogy(SNR,vector1,'o-','Color', colorstring(1),'LineWidth',1,'MarkerSize',4,'MarkerFaceColor','b')
hold on;grid on; ylim([10^(-5) 0.446]);
semilogy(SNR.*b,vector2,'s-','Color', colorstring(2),'LineWidth',1,'MarkerSize',4,'MarkerFaceColor','m')%7
semilogy(SNR.*c,vector3,'h-','Color', colorstring(3),'LineWidth',1,'MarkerSize',4,'MarkerFaceColor','r')%11
legend('User 1 \alpha_1 = 0.70','User 2 \alpha_2 = 0.20','User 3 \alpha_3 = 0.05','Location','southwest');
title('BER for mth user in Rayleigh fading channel');
xlabel('SNR');
ylabel('BER');
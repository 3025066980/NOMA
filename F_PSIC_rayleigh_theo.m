clc; clear variables;
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
vector1=zeros([1 length(SNR)]);
vector2=zeros([1 length(SNR)]);
vector3=zeros([1 length(SNR)]);
vector4=zeros([1 length(SNR)]);
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
               omega=sqrt( 1+( (2*(k+1)) / (beta^2*gamma(h)*1) ) );
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
               Y=mtx2(i,:);
               beta=sum(Y)/sqrt(P);
               omega=sqrt( 1+( (2*(k+1)) / (beta^2*gamma(h)*1) ) );
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
               Y=mtx3(i,:);
               beta=sum(Y)/sqrt(P);
               omega=sqrt( 1+( (2*(k+1)) / (beta^2*gamma(h)*1) ) );
               p1=(ita*(omega-1)*(-1)^(t+k)*nchoosek((M-m),t)*nchoosek((m+t-1),k)) / (2^(M-m+1)*(k+1)*omega);
               sum1_1=sum1_1+p1;
           end
           sum1_2=sum1_2+sum1_1;
        end
        sum1_3=sum1_3+sum1_2;
    end   
    vector3(h)=sum1_3;
end
for h=1:length(SNR)
    sum1_3=0;
    m=4;
    for i=1:ls4
        sum1_2=0;
        for t=0:(M-m)
           sum1_1=0;
           for k=0:(m+t-1)
               ita=factorial(M)/(factorial(M-m)*factorial(m-1));
               Y=mtx4(i,:);
               beta=sum(Y)/sqrt(P);
               omega=sqrt( 1+( (2*(k+1)) / (beta^2*gamma(h)*1) ) );
               p1=(ita*(omega-1)*(-1)^(t+k)*nchoosek((M-m),t)*nchoosek((m+t-1),k)) / (2^(M-m+1)*(k+1)*omega);
               sum1_1=sum1_1+p1;
           end
           sum1_2=sum1_2+sum1_1;
        end
        sum1_3=sum1_3+sum1_2;
    end   
    vector4(h)=sum1_3;
end
colorstring = 'bmrgk';
figure(2)
semilogy(SNR,vector1,'o-','Color', colorstring(1),'LineWidth',1,'MarkerSize',4,'MarkerFaceColor','b')
hold on;grid on; ylim([10^(-5) 1]);
semilogy(SNR,vector2,'s-','Color', colorstring(2),'LineWidth',1,'MarkerSize',4,'MarkerFaceColor','m')
semilogy(SNR,vector3,'h-','Color', colorstring(3),'LineWidth',1,'MarkerSize',4,'MarkerFaceColor','r')
semilogy(SNR,vector4,'d-','Color', colorstring(4),'LineWidth',1,'MarkerSize',4,'MarkerFaceColor','g')
legend('User 1 \alpha_1 = 0.75','User 2 \alpha_2 = 0.18','User 3 \alpha_3 = 0.05','User 4 \alpha_4 = 0.02','Location','southwest');
title('BER for mth user in Rayleigh fading channel');
xlabel('SNR');
ylabel('BER');
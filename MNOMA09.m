clc
clear 
M = 20;     %- Numero total de usuarios
aux_a = (2)*((1/4).^(0:M-1));
a = aux_a/sum(aux_a); %- a Normalizado
P = 1;
alpha = sqrt(a*P);
SNRdB = 0:50;
SNR = 10.^(SNRdB/10);
b=zeros(1,51);
BER=zeros(51,20);
%% Simulacion
NS = 1e5;   %- Numero de simulaciones
xi = rand(NS,M)>0.5;
XT = 2*xi-1;
XT = alpha.*XT;
YT = sum(XT,2);
PNoise = (P)./SNR;
media=zeros(length(SNR),M);
% i=15;
for i = 1:length(SNR)%SNR
    XHAT = zeros(NS,1);
    for k = 1:M%ususarios
        N = wgn(NS,1,10*log10(PNoise(i)));
        RT = YT+N-sum(XHAT,2);
        xhat = RT>0;
        XHAT(:,k) = alpha(k)*(2*xhat-1);
        media(i,k)=mean(XHAT(:,k));
        EM(i).N(:,k) = sum(XHAT(:,1:k) - XT(:,1:k),2);
        BER(i,k) = biterr(xi(:,k),xhat)/NS;
    end
end
%%
mean(media(20,:))
% stem(media(1,12))
%%
clc;% close all;
% S =19;%SNRb
I =8;%User
for S=1:length(SNR)
    % histogram(EM(N).N(:,I))
    b(S)=var(EM(S).N(:,I));%% simulated VAR
end
% g=a(1);
% for k=2:I
%     g=g-a(k);
% end
p=0.3109+0.4512./(1:M);
f=p(I)*exp(-0.18.*SNRdB);%user 11
% f=0.75*exp(-0.115.*SNRdB);%user 11
% f=0.6*exp(-0.115.*SNRdB);%user 22
% f=0.50*exp(-0.14.*SNRdB);%user 33
% f=0.50*exp(-0.18.*SNRdB);%user 44
% f=0.42*exp(-0.18.*SNRdB);%user 55
% f=0.39*exp(-0.19.*SNRdB);%user 66

% f=exp(-0.271-0.114.*SNRdB);%user 1
% f=exp(-0.558-0.115.*SNRdB);%user 2
% f=exp(-0.965-0.118.*SNRdB);%user 3
% f=exp(-1.260-0.127.*SNRdB);%user 4
% f=exp(-1.401-0.141.*SNRdB);%user 5
% f=exp(-1.397-0.157.*SNRdB);%user 6
% hold off
plot(SNRdB,f),hold on;
stem(b),hold off;
xlabel('SNR')
title(['VAr for ',num2str(I),'th user'])
ylabel('VAR')
legend('VAR')
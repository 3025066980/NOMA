clc
clear 
M = 20;%- Numero total de usuarios
aux_a = (2)*((1/4).^(0:M-1));
a = aux_a/sum(aux_a); %- a Normalizado
P = 1;
alpha = sqrt(a*P);
SNRdB = 0:50;
SNR = 10.^(SNRdB/10);
b=zeros(1,51);
BER=zeros(51,20);
%% Simulacion
NS = 1e5;%- Numero de simulaciones
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
clc;
% for I=7:7
I =7;%User
for S=1:length(SNR)
    % histogram(EM(N).N(:,I))
    b(S)=var(EM(S).N(:,I));%% simulated VAR
end
stem(b),hold on;
gg=sum(a(I+1:M))/a(I);%coeficiente a
hh= 0.1437*exp(0.5269*log(alpha/1.732))-0.2012;%coeficiente b
% end
f=gg*exp(hh(I).*SNRdB);%Funcion varianza
plot(f)%,hold on;
xlabel('SNR')
title(['Var for ',num2str(I),'th user'])
ylabel('VAR')
legend('VAR')
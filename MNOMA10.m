clc
clear 
M = 20;%- Numero total de usuarios
aux_a = (2)*((1/5).^(0:M-1));
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
% for I=1:20
I =1;%User
for S=1:length(SNR)
    % histogram(EM(N).N(:,I))
    b(S)=var(EM(S).N(:,I));%% simulated VAR
end
stem(SNRdB,b),hold on;
% end
gg=0.383.*alpha+0.3706;
hh= (0.1489.*a + 0.002821) ./ (a + 0.01485);
f=gg(I)*exp(-hh(I).*SNRdB);%Funcion varianza
plot(SNRdB,f)%,hold on;
xlabel('SNR')
title(['Var for ',num2str(I),'th user'])
ylabel('VAR')
legend('VAR')
% gg=[0.7046 0.5432 0.4395 0.3999 0.3833 0.376 0.3727 0.3712 0.3705 0.3702 0.3702 0.37 0.37 0.37 0.37 0.37 0.37 0.37 0.37 0.37];
% hh=[0.157 0.1436 0.158 0.1744 0.1831 0.187 0.1888 0.1895 0.1898 0.19 0.19 0.1901 0.1901 0.1901 0.1901 0.1901 0.1901 0.1901 0.1901 0.1901];

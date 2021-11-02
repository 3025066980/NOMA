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
I =20;%User
for S=1:length(SNR)
    % histogram(EM(N).N(:,I))
    b(S)=var(EM(S).N(:,I));%% simulated VAR
end
stem(b),hold on;
% end
% gg=[0.6727 0.5240 0.4186 0.3736 0.3517 0.3405 0.3348 0.3320 0.3306 0.3298 0.3295 0.3293 0.3292 0.3292 0.3292 0.3291 0.3291 0.3291 0.3291 0.3291];
% hh=[-0.1123 -0.1173 -0.1429 -0.1696 -0.1847 -0.1921 -0.1956 -0.1974 -0.1982 -0.1986 -0.1988 -0.1989 -0.1990 -0.1990 -0.1990 -0.1990 -0.1990 -0.1990 -0.1990 -0.1990];
% hh=[0.1123 0.1173 0.1429 0.1696 0.1847 0.1921 0.1956 0.1974 0.1982 0.1986 0.1988 0.1989 0.1990 0.1990 0.1990 0.1990 0.1990 0.1990 0.1990 0.1990];
gg=0.4066.*alpha+0.3298;
hh= (0.109.*a + 0.00491) ./ (a + 0.02478);
f=gg(I)*exp(-hh(I).*SNRdB);%Funcion varianza
plot(f)%,hold on;
xlabel('SNR')
title(['Var for ',num2str(I),'th user'])
ylabel('VAR')
legend('VAR')
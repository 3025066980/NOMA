clc
clear 
M = 30;%- Numero total de usuarios
aux_a = (2)*((1/4).^(0:M-1));
a = aux_a/sum(aux_a); %- a Normalizado
P = 1;
alpha = sqrt(a*P);
SNRdB = 0:50;
SNR = 10.^(SNRdB/10);
b=zeros(1,51);
BER=zeros(51,20);
%% Simulacion
NS = 1e4;%- Numero de simulaciones
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
% u=zeros(1,M);
% for I=1:M
I =20;%User
for S=1:length(SNRdB)
    % histogram(EM(N).N(:,I))
    b(S)=var(EM(S).N(:,I));%% simulated VAR
end
% u(I)=b(1);
stem(SNRdB,b),hold on;
% end
gg=0.4065.*alpha+0.3287;
hh= (0.1084.*a + 0.005051) ./ (a + 0.02556);
f=gg(I)*exp(-hh(I).*SNRdB);%Funcion varianza
plot(SNRdB,f)%,hold on;
xlabel('SNR')
title(['Var for ',num2str(I),'th user'])
ylabel('VAR')
legend('VAR')
%%
rr=zeros(1,M);
for k=1:M
rr(k)=sum(alpha(1:k));
end
stem(-rr*a(1))
%%
N=15;%SNR
I=8;%usuario
histogram(EM(N).N(:,I))
% gg=[0.6716 0.5225 0.4171 0.372 0.3529 0.3417 0.3331 0.3302 0.3288 0.3281 0.3277 0.3277 0.3275 0.3274 0.3274 0.3274 0.3274 0.3274 0.3274 0.3274];
% hh=[0.1118 0.117 0.1427 0.1691 0.1854 0.1927 0.195 0.1967 0.1975 0.1979 0.1981 0.1981 0.1983 0.1983 0.1983 0.1983 0.1983 0.1983 0.1983 0.1983];
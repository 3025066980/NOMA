%Imperfect SIC
clc
clear 
%close all
M =15;     %- Numero total de usuarios
aux_a = (2)*((1/4).^(0:M-1));%coeficientes de canal
a = aux_a/sum(aux_a); %- coeficienes Normalizados
P = 1;
alpha = sqrt(a*P);
alphasum=0:40;
SNRdB = 0:50;
SNR = 10.^(SNRdB/10);
NS = 1e5;   %- Numero de simulaciones
%% Simulacion   %- Numero de simulaciones
xi = rand(NS,M)>0.5;
XTt = 2*xi-1;
XT = alpha.*XTt;
% xxhat=zeros(length(SNR),M);
YT = sum(XT,2);
PNoise = (P)./SNR;
for i = 1:length(SNR)
    XHAT = zeros(NS,1);
    for k = 1:M
        N = wgn(NS,1,10*log10(PNoise(i)));
        RT = YT+N-sum(XHAT,2);
        xhat = RT>0;
%         xxhat(:,k)=2*xhat-1;%alt
        XHAT(:,k) = alpha(k)*(2*xhat-1);
        EM(i).N(:,k) = sum(XHAT(:,1:k) - XT(:,1:k),2);
        BER(i,k) = biterr(xi(:,k),xhat)/NS;
    end
%     media(i,:)=mean(xxhat);%alt
end
%% Probabilidades de error 
figure()
semilogy(SNRdB,BER)
%% Histogramas
clc
% for I=1:M%para recorrer los usuarios
a=0;
N = 11;%SNR
I = 13;%User
% histogram(EM(N).N(:,I),100)
var(EM(N).N(:,I))%% VAR theo

al=0;
for l=1:I-1
   al=(alpha(l)^2)*BER(N,l)+al;
end
4*al
disp('------------------')
%% Means
%mean of x
mean(XTt)
%mean of x hat
mean(media)
%% 2nd term
clc;
NS = 1e5;   %- Numero de simulaciones
xi = rand(NS,M)>0.5;
XTt = 2*xi-1;
XT = alpha.*XTt;
% xxhat=zeros(length(SNR),M);
YT = sum(XT,2);
PNoise = (P)./SNR;
i=15;%snr position
m=10;%m-th user
XHAT = zeros(NS,1);
for k = 1:M
    N = wgn(NS,1,10*log10(PNoise(i)));
    RT = YT+N-sum(XHAT,2);
    xhat = RT>0;
    xxhat(:,k)=2*xhat-1;
    XHAT(:,k) = alpha(k)*(2*xhat-1);
    EM(i).N(:,k) = sum(XHAT(:,1:k) - XT(:,1:k),2);
    BER(i,k) = biterr(xi(:,k),xhat)/NS;
end
all=0;
for t=1:m-1
    for j=1:m-1
        al=alpha(t)*alpha(j).*xxhat(:,t).*xxhat(:,j);
        all=all+al;
    end
end
g=mean(all);
alfa=0;
for e=1:m-1
   alfa=alpha(e)^2+alfa;
end
fprintf('2nd term is %d \nsum alfa is %d\n',g,alfa)
%% 3rd term
clc;
NS = 1e5;   %- Numero de simulaciones
xi = rand(NS,M)>0.5;
XTt = 2*xi-1;
XT = alpha.*XTt;
% xxhat=zeros(length(SNR),M);
YT = sum(XT,2);
PNoise = (P)./SNR;
i=45;%snr position
m=6;%m-th user
XHAT = zeros(NS,1);
for k = 1:M
    N = wgn(NS,1,10*log10(PNoise(i)));
    RT = YT+N-sum(XHAT,2);
    xhat = RT>0;
    xxhat(:,k)=2*xhat-1;
    XHAT(:,k) = alpha(k)*(2*xhat-1);
    EM(i).N(:,k) = sum(XHAT(:,1:k) - XT(:,1:k),2);
    BER(i,k) = biterr(xi(:,k),xhat)/NS;
end
all=0;
for t=1:m-1
    for j=1:m-1
        al=alpha(t)*alpha(j).*XTt(:,t).*xxhat(:,j);
        all=all+al;
    end
end
g=mean(all);
alfa=0;
for e=1:m-1
   alfa=alpha(e)^2+alfa;
end
fprintf('3nd term is %d \nsum alfa is %d\n',g,alfa)
%% 4rd term
clc;
NS = 1e5;   %- Numero de simulaciones
xi = rand(NS,M)>0.5;
XTt = 2*xi-1;
XT = alpha.*XTt;
% xxhat=zeros(length(SNR),M);
YT = sum(XT,2);
PNoise = (P)./SNR;
i=30;%snr position
m=6;%m-th user
XHAT = zeros(NS,1);
for k = 1:M
    N = wgn(NS,1,10*log10(PNoise(i)));
    RT = YT+N-sum(XHAT,2);
    xhat = RT>0;
    xxhat(:,k)=2*xhat-1;
    XHAT(:,k) = alpha(k)*(2*xhat-1);
    EM(i).N(:,k) = sum(XHAT(:,1:k) - XT(:,1:k),2);
    BER(i,k) = biterr(xi(:,k),xhat)/NS;
end
all=0;
for t=1:m-1
    for j=1:m-1
        al=alpha(t)*alpha(j).*xxhat(:,j).*XTt(:,t);
        all=all+al;
    end
end
g=mean(all);
alfa=0;
for e=1:m-1
   alfa=alpha(e)^2+alfa;
end
fprintf('4th term is %d \nsum alfa is %d\n',g,alfa)

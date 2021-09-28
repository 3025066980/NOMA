%perfect SIC
clc
clear 
%close all

M = 10; %- Numero total de usuarios
m = 5; %- Usuario que estoy codificando

aux_a = (3/4)*((1/4).^(0:M-1));
a = aux_a/sum(aux_a); %- a Normalizado
P = 1;
alpha = sqrt(a*P);

% Distribucion de I
NS = 1e4;
S = 2*(rand(NS,M)>0.5)-1;
S_m = S(:,m+1:end);
I_sim = S_m*alpha(m+1:end)';
sum(alpha(m+1:end))
figure(1)
histogram(I_sim,500)
%% Calculo segun el articulo
J = 2^(M-m);
A = 0:J-1;
B = ones(J,1);

alpha_cut = alpha(m:end);
lambda = [B  2*de2bi(A')-1];
lambdaj = lambda*alpha_cut';

SNRdB = 0:40;
SNR = 10.^(SNRdB/10);
for i = 1:length(SNR)
    Pe(i)= 1/J * sum(qfunc(lambdaj*sqrt(SNR(i))/sqrt(P)));
end
figure(1)
semilogy(SNRdB,Pe)
ylim([1e-5 1])



%% Calculo si consideramos Gaussiana
varI = 1*sum(alpha(m+1:end).^2);

for i = 1:length(SNR)
    PeT(i) = qfunc(sqrt(alpha(m)^2/(varI+(P)/SNR(i))));
end

figure(1)
hold on
semilogy(SNRdB,PeT,'*')

%% Simulacion
NS = 1e4; %- Numero de simulaciones
xi = rand(NS,M)>0.5;
XT = 2*xi-1;
XT = alpha.*XT;
XT_m = XT(:,m:end);

YT = sum(XT_m,2);

PNoise = (P)./SNR;
for i = 1:length(SNR)
    N = wgn(NS,1,10*log10(PNoise(i)));
    RT = YT+N;
    xhat = RT>0;
    BER(i) = biterr(xi(:,m),xhat)/NS;
end

figure(1)
semilogy(SNRdB,BER,'--')




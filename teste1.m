clc
clear 
%close all

M = 5; %- Numero total de usuarios
m = 2; %- Usuario que estoy codificando

aux_a = (1)*((1/4).^(0:M-1));
a = aux_a/sum(aux_a); %- a Normalizado
P = 1;
alpha = sqrt(a*P);

alim = sum(alpha(m+1:end));
beta = linspace(-alim,alim,1e4);

SNRdB = 0:40;
SNR = 10.^(SNRdB/10);

for i = 1:length(SNR)
    cmas = (alpha(m) + alim)/sqrt(P/SNR(i));
    cmenos = (alpha(m) - alim)/sqrt(P/SNR(i));    
    PeTA(i) = sqrt(P/SNR(i))/(2*alim)*(cmas * qfunc(cmas) - cmenos * qfunc(cmenos) - (exp(-1*(cmas^2)/2) - exp(-1*(cmenos^2)/2))/sqrt(2*pi));
end

figure(1)
%semilogy(SNRdB,PeT,'*')
%hold on
semilogy(SNRdB,PeTA,'--')

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
hold on
semilogy(SNRdB,Pe)
ylim([1e-5 1])
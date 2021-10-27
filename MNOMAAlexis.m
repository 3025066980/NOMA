clc
clear 
M = 20;
aux_a = (2)*((1/4).^(0:M-1));
a = aux_a/sum(aux_a);
P = 1;
alpha = sqrt(a*P);
SNRdB = 0:50;
SNR = 10.^(SNRdB/10);
% Función Q, P(e)
lamda =SNR;
I = 10;
deviation = sqrt(P./lamda);
al = 0;
for i=I+1:M
    al = al + alpha(i);
end
b = (alpha(I) + al)./deviation;
c = (alpha(I) - al)./deviation;
P_e = b.*qfunc(b) - c.*qfunc(c) - (exp((-b.^2)/2) + exp((-c.^2)/2))/sqrt(2*pi);
plot(P_e)
xlabel("SNR")
ylabel("P(e)")
title(['P(e) for ',num2str(I),'th user'])
grid on
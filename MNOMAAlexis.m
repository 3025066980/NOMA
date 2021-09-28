clc
clear 
%close all

M = 20;     %- Numero total de usuarios
aux_a = (2)*((1/4).^(0:M-1));
a = aux_a/sum(aux_a); %- a Normalizado
P = 1;
alpha = sqrt(a*P);

SNRdB = 0:50;
SNR = 10.^(SNRdB/10);
%% Simulacion
NS = 1e6;   %- Numero de simulaciones
xi = rand(NS,M)>0.5;
XT = 2*xi-1;
XT = alpha.*XT;

YT = sum(XT,2);
PNoise = (P)./SNR;
for i = 1:length(SNR)
    XHAT = zeros(NS,1);
    for k = 1:M
        N = wgn(NS,1,10*log10(PNoise(i)));
        RT = YT+N-sum(XHAT,2);
        xhat = RT>0;
        XHAT(:,k) = alpha(k)*(2*xhat-1);
        EM(i).N(:,k) = sum(XHAT(:,1:k) - XT(:,1:k),2);
        BER(i,k) = biterr(xi(:,k),xhat)/NS;
    end
end
%% Probabilidades de error 
figure()
semilogy(SNRdB,BER)

%% Histrogramas
I = 20; %N°usuario

AM = XT-XHAT;
AM(find(AM)) = 1;
AM1 = sum(AM,1);
P_e = AM1/NS;

%sum
T_var = 4*alpha.*P_e;
varA = 0;
P_eMm = zeros([1 20]);
P_emM = zeros([1 20]);
P_e_1 = 0;
P_e_0 = 0;
for u = 1:I
    varA(u+1) = sum(T_var(1:u)); 
    for v = 1:1000000
        %S=1 y -S=-1
        if XT(v,u)>0 && XHAT(v,u)<0
            P_e_1 = P_e_1 + 1;
        end 
        %S=-1 y -S=1
        if XT(v,u)<0 && XHAT(v,u)>0
            P_e_0 = P_e_0 + 1;
        end
    end
    P_eMm(u) = P_e_1; 
    P_emM(u) = P_e_0; 
    P_e_1=0;
    P_e_0=0;
end
P1_1 = P_eMm/NS
P_11 = P_emM/NS
%varA
%%
for usuario=1:M
    usuario;
    N;
    T_var = var(EM(N).N(:,usuario)); %teoric
    
    DT = 0;
    real = XT./alpha;
    real2 = XHAT./alpha;
    for a=1:NS
        if real2(a,usuario) ~= real(a,usuario)
            DT = DT + 1;
        end
    end
    s_var = 0;
    for b = 1:usuario-1
        s_var = alpha(b)*(DT/NS);
    end 
    S_var = 4*(s_var);
    %almacenamos
    comp(usuario,1) = T_var;
    comp(usuario,2) = S_var;
end
comp
%% Histogramas
clc
for I=1:M%para recorrer los usuarios
N = 1;%SNR
% I = 1;%User
histogram(EM(N).N(:,I),100)
var(EM(N).N(:,I))%% VAR theo
al=0;
for l=1:I-1
   al=alpha(l)^2*BER(N,l)+al;
end
4*al
disp('------------------')
end

%%
tabla = zeros(20,50);
for I=1:20
    for N=1:50
        I;
        N;
        tabla(I,N) = var(EM(N).N(:,I));
    end
end
%xlswrite('SNR & var.xlsx', tabla,'Hoja1','C3');
%%
i = 1;%numero de usuario
a = 0.00007091*i^4 - 0.00355893*i^3 + 0.06106485*i^2 - 0.39515201*i + 1.11382883;
b = -0.00000950*i^4 + 0.00041726*i^3 - 0.00548675*i^2 + 0.01347378*i - 0.12148400;
func = a*exp(b.*SNRdB);
figure()
grid on
hold on
plot(SNRdB,func,':blue')
plot(SNRdB(1:50),tabla(i,1:50),'red')%% VAR theo
legend('Simulado','Real')
title(['Usuario', i])

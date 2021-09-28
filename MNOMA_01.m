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
N = 50; % SNR
I = 20; %N°usuario
%figure()
%histogram(EM(N).N(:,I),100)


%%Comparamos datos (xi y xhat)
comp = zeros(M,2);
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
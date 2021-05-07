clc; clear variables; 
% close all;

df = 800; dn = 200; 	%Distances
% df = 800; dn = 500;
eta = 4;	%Path loss exponent
N = 10^5;

%Rayleigh fading coefficients
hf = sqrt(df^-eta)*(randn(N,1)+1i*randn(N,1))/sqrt(2);
hn = sqrt(dn^-eta)*(randn(N,1)+1i*randn(N,1))/sqrt(2);

%Channel gains
gf = (abs(hf)).^2;
gn = (abs(hn)).^2;

%Transmit power
Pt = -60:5:60;	%in dBm
pt = (10^-3)*db2pow(Pt); %in watts

BW = 10^6;	%bandwidth
%Noise powers
No = -174 + 10*log10(BW);	%in dBm
no = (10^-3)*db2pow(No);	%in watts

pf = zeros(1,length(pt));
pn = zeros(1,length(pt));
Rf = zeros(1,length(pt));
Rn = zeros(1,length(pt));

rf = 0.5; rn = 0.5;	%Target rates

for u = 1:length(pt)
    Cf = log2(1 + gf*pt(u)/no);	%Rate of far user
    Cn = log2(1 + gn*pt(u)./(gf*pt(u)+no));	%Rate of near user
    
    for k = 1:N
        if (Cf(k) < rf) || (Cn(k) < rn) %Outage condition for far user
            pf(u) = pf(u) + 1;
        end
        if Cn(k) < rn	%Outage condition for near user
            pn(u) = pn(u)+1;
        end
        
    end
end

%figure;
semilogy(Pt,pf/N,'-*r','linewidth',1.5); hold on; grid on;
semilogy(Pt,pn/N,'-*b','linewidth',1.5);

xlabel('Transmit SNR (dB)');
ylabel('Outage probability');
title('Outage probability of uplink NOMA');
legend('Far user','Near user');


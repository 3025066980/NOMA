clc; clear variables;
N = 10^5;
SNR = 0:40;
snr = 10.^(SNR/10);
P=1;
a1 = 0.75; a2 = 0.18; a3=0.05; a4=0.02;
ber1=zeros(1,length(SNR));
ber2=zeros(1,length(SNR));
ber3=zeros(1,length(SNR));
ber4=zeros(1,length(SNR));
a=randn(1);
for u = 1:length(SNR)
    %Generate random binary data
    data1 = randi([0 1],1,N);  %Data bits of user 1
    data2 = randi([0 1],1,N);  %Data bits of user 2
    data3 = randi([0 1],1,N);  %Data bits of user 3
    data4 = randi([0 1],1,N);  %Data bits of user 4
    %Do BPSK modulation of data
    xmod1 = 2*data1 - 1;
    xmod2 = 2*data2 - 1;
    xmod3 = 2*data3 - 1;
    xmod4 = 2*data4 - 1;
    %Rayleigh coefficent
    h1=zeros(2,N);
    h1(1,:)=raylrnd(1/sqrt(8),1,N);
    h1(2,:)=h1(1,:).^2;
    h1=h1';
    H1 = sortrows(h1,2);
    H1=H1';
    
    h2=zeros(2,N);
    h2(1,:)=raylrnd(1/sqrt(4),1,N);
    h2(2,:)=h2(1,:).^2;
    h2=h2';
    H2 = sortrows(h2,2);
    H2=H2';
    
    h3=zeros(2,N);
    h3(1,:)=raylrnd(1/sqrt(2),1,N);
    h3(2,:)=h3(1,:).^2;
    h3=h3';
    H3 = sortrows(h3,2);
    H3=H3';
    
    h4=zeros(2,N);
    h4(1,:)=raylrnd(1/sqrt(1),1,N);
    h4(2,:)=h4(1,:).^2;
    h4=h4';
    H4 = sortrows(h4,2);
    H4=H4';
    
    x = sqrt(P)*(sqrt(a1)*xmod1 + sqrt(a2)*xmod2 + sqrt(a3)*xmod3 + sqrt(a4)*xmod4);
    N0 = 1/10^(-SNR(u)/10);
    %Received signals
    y1 = H1(1,:).*x + 1/sqrt(N0)*(randn(1,N));
    y2 = H2(1,:).*x + 1/sqrt(N0)*(randn(1,N));
    y3 = H3(1,:).*x + 1/sqrt(N0)*(randn(1,N));
    y4 = H4(1,:).*x + 1/sqrt(N0)*(randn(1,N));
    %Equalize
    eq1 = y1./H1(1,:);
    eq2 = y2./H2(1,:);
    eq3 = y3./H3(1,:);
    eq4 = y4./H4(1,:);
    %AT USER 1------------------------
    %Direct decoding of x from y1
    x1_hat = ones(1,N);         
    x1_hat(eq1 < 0) = 0;         %Final bits for user 1
    %AT USER 2------------------------
    rem2=eq2-sqrt(a1*P)*xmod1;
    %Decode x2 from rem without signal x1
    x2_hat = zeros(1,N);
    x2_hat(rem2>0) = 1;          %Final bits for user 2
    %AT USER 3-------------------------
    rem3 = eq3 - sqrt(a1*P)*xmod1 - sqrt(a2*P)*xmod2;
    %Decode x3 from rem
    x3_hat = zeros(1,N);
    x3_hat(rem3>0) = 1;          %Final bits for user 3
    %AT USER 4------------------------
    rem4=eq4 - sqrt(a1*P)*xmod1 - sqrt(a2*P)*xmod2 - sqrt(a3*P)*xmod3;
    x4_hat = zeros(1,N);
    x4_hat(rem4>0) = 1;
       
    ber1(u) = biterr(data1,x1_hat)/N;
    ber2(u) = biterr(data2,x2_hat)/N;
    ber3(u) = biterr(data3,x3_hat)/N;
    ber4(u) = biterr(data4,x4_hat)/N;
end
colorstring = 'bmrgk';
figure(2)
semilogy(SNR, ber1,'+--','Color', colorstring(1), 'linewidth', 1); 
ylim([10^(-5) 1])
hold on; grid on;
semilogy(SNR, ber2,'+--','Color', colorstring(2), 'linewidth', 1);
semilogy(SNR, ber3,'+--','Color', colorstring(3), 'linewidth', 1); 
semilogy(SNR, ber4,'+--','Color', colorstring(4), 'linewidth', 1); 
xlabel('SNR');
ylabel('BER');
clc; clear variables;
user=input('Number of users: ');
histograma=input('Plot Histogram? ');
graph=input('Plot graphs? ');
reqsnr=input('Histogram in SNR= ');
a1=3/4; N = 1e5; SNR = 0:40; P=1;
powcoef=zeros(1,user);
powcoef(1)=a1;
for i=2:user-1
    powcoef(i)=a1*(1-sum(powcoef(1:i-1)));
end
powcoef(user)=1-sum(powcoef);
% powcoef=[0.70,0.25 0.05]; N=1e5; SNR=0:40;P=1;
ber=zeros(user,N);
biter=zeros(user, length(SNR));
ems=cell(length(SNR),1);
for z=1:length(SNR)
    ems{z}=zeros(user, N);
end
for u = 1:length(SNR)
    x=zeros(user,N);
    xmod=zeros(user,N);
    xx=zeros(1,N);
    noise=zeros(user,N);
    y=zeros(user,N);
    for i=1:user
        x(i,:) = randi([0 1],1,N);%Create 0 and 1
    end
    for j=1:user
        xmod(j,:) = 2*x(j,:)-1;%Bpsk bits
    end
    for k=1:user
        xx=xx+sqrt(powcoef(k)*P)*xmod(k,:);%Superposition
    end
    for l=1:user%creating noise for each user
        f=0;
        for z=l:user
            f=f+powcoef(z);
        end
        noise(l,:) = wgn(1,length(xx),10*log10(f/(10^(SNR(u)/10))));
    end
    for m=1:user%Applying user's noise to superpositioned signal
       y(m,:)=xx+noise(m,:); 
    end
    for n=1:user
        if n==1
            x1_hat = ones(1,N);
            x1_hat(y(1,:) < 0) = 0;
            ber(1,:)=x1_hat;
        else
            for o=1:n-1
                x_hat = ones(1,N);%SIC process
                x_hat(y(n,:) < 0) = -1;%Estiamate
                y(n,:)=y(n,:)-sqrt(powcoef(o)*P)*x_hat;%Subtracting m-1 user to m
                ems{u}(n,:)=ems{u}(n,:)+(xmod(n,:)-x_hat);
            end
            x_zero=ones(1,N);%BPSK demodulation
            x_zero(y(n,:)<0)=0;
            ber(n,:)=x_zero;
        end
    end
    for p=1:user
       biter(p,u)=biterr(x(p,:),ber(p,:))/N;
    end
end
colorstring = 'bmryk';
if graph==1
    figure(1)
    for v=1:user
        txt = ['U_{',num2str(v),'} \alpha=',num2str(powcoef(v))];
        semilogy(SNR, biter(v,:),'+--', 'linewidth', 1,'DisplayName',txt);
        grid on; hold on;
    end
    title('BER in awgn channel')
    xlabel('SNR(in dB)')
    ylabel('BER')
    legend show
end
if histograma==1
    figure(2)
    title(['e_m for SNR=',num2str(reqsnr)])
    for b=2:user
        txt = ['e_{',num2str(b),'}'];
        histogram(ems{reqsnr}(b,:),'DisplayName',txt)
        legend show
        title(['e_m for SNR=',num2str(reqsnr)])
        pause
        hold on;
    end
end
close all;
clear;

fc = 800*10^6;    %carrier frequency
fb = 2*10^6;      %bit frequency is twice of the symbol frequency
Tb = 1/fb;        %bit time period
Th = -135;        %given value of outage threshold
Ts = 1/(4*fc);    %sampling time should be less than the twice of the carrier time period
D = 0.1;          %given distance between tx and rx in km
j = 1;            %iteration number

total_syms = 1000;           %total no of symbols to be taken
total_bits = 2*total_syms;     %total no of bits
data1 = round(rand(1,total_bits));  %generating random bits of 0 and 1
data2 = 2*data1 - 1;                %converting 0 to -1 and 1 to 1
data3 = reshape(data2,2,length(data1)/2); %here we separate the inphase and quadphase bits
t = 0:Ts:Tb-Ts;  
tl = length(t);
[tx_bits, Xi, Xq] = modulate(t,fc,data3);    %modulating the transmit bits 
L = size(tx_bits,2);
b = fir1(50,0.05,'low');    %here we create low pass filter for demodulation

for A = 10^(-3):10^(-3):10^(-2)
    tx_signal = A*tx_bits;             %amplitude of transmitting signal
    %here we create channel gains for three different channels
    cg_chnl1 = sqrt(0.5)*randn(1,total_syms) + sqrt(0.5)*1i*randn(1,total_syms); %here we generate complex gaussian channel based on rayleigh parameter
    h1 = repelem(cg_chnl1,tl);    %here we create slow fading channel with same channel gain for one bit period
    cg_chnl2 = sqrt(0.5)*randn(1,total_syms) + sqrt(0.5)*1i*randn(1,total_syms);   %here we generate complex gaussian channel based on rayleigh parameter
    h2 = repelem(cg_chnl2,tl);     %here we create slow fading channel with same channel gain for one bit period
    cg_chnl3 = sqrt(0.5)*randn(1,total_syms) + sqrt(0.5)*1i*randn(1,total_syms);    %here we generate complex gaussian channel based on rayleigh parameter
    h3 = repelem(cg_chnl3,tl);     %here we create slow fading channel with same channel gain for one bit period
    
    pathloss_dB = 128.1 + 37.6*log10(D);
    pathloss = 10^-(pathloss_dB/20);        %pathloss

    
    P = 10^((-100-30)/10);          %calculation of power
    sigma = sqrt(P/2);
    %here we randomly generate the noise for three channels
    noise1 = sigma*randn(1,L)+1i*sigma*randn(1,L);
    noise2 = sigma*randn(1,L)+1i*sigma*randn(1,L);
    noise3 = sigma*randn(1,L)+1i*sigma*randn(1,L);
    
    signal1 = pathloss*h1.*tx_signal;      %signal without noise
    rx_signal1 = signal1 + noise1;         %signal with noise
    signal2 = pathloss*h2.*tx_signal;
    rx_signal2 = signal2 + noise2;
    signal3 = pathloss*h3.*tx_signal;
    rx_signal3 = signal3 + noise3;

    conj_h1 = conj(h1);
    Yh1 = real(conj_h1.*rx_signal1./abs(h1));
    conj_h2 = conj(h2);
    Yh2 = real(conj_h2.*rx_signal2./abs(h2));
    conj_h3 = conj(h3);
    Yh3 = real(conj_h3.*rx_signal3./abs(h3));

    inphase_rx = [];
    quadphase_rx = [];
    outage = 0;
    for i = 1 : total_syms
        pow_sym1_dBm = pow2db((1/tl)*sum(real(rx_signal1((i-1)*tl+1:i*tl)).^2))+30;
        pow_sym2_dBm = pow2db((1/tl)*sum(real(rx_signal2((i-1)*tl+1:i*tl)).^2))+30;
        pow_sym3_dBm = pow2db((1/tl)*sum(real(rx_signal3((i-1)*tl+1:i*tl)).^2))+30;
        [~, ind] = max([pow_sym1_dBm,pow_sym2_dBm,pow_sym3_dBm]);
        %Keeping only the maximum power signals and demodulating them
        %also computing power for outage probability
        if(ind==1)
            sym_power = (1/tl)*sum(real(signal1((i-1)*tl+1:i*tl)).^2);
            sym_power_dBm = pow2db(sym_power)+30;
            inphase_rx((i-1)*tl+1:i*tl) = Yh1((i-1)*tl+1:i*tl).*Xi;
            quadphase_rx((i-1)*tl+1:i*tl) = Yh1((i-1)*tl+1:i*tl).*Xq;
        elseif(ind==2)
            sym_power = (1/tl)*sum(real(signal1((i-1)*tl+1:i*tl)).^2);
            sym_power_dBm = pow2db(sym_power)+30;
            inphase_rx((i-1)*tl+1:i*tl) = Yh2((i-1)*tl+1:i*tl).*Xi;
            quadphase_rx((i-1)*tl+1:i*tl) = Yh2((i-1)*tl+1:i*tl).*Xq;
        else
            sym_power = (1/tl)*sum(real(signal1((i-1)*tl+1:i*tl)).^2);
            sym_power_dBm = pow2db(sym_power)+30;
            inphase_rx((i-1)*tl+1:i*tl) = Yh3((i-1)*tl+1:i*tl).*Xi;
            quadphase_rx((i-1)*tl+1:i*tl) = Yh3((i-1)*tl+1:i*tl).*Xq;
        end
        if(sym_power_dBm < Th)
            outage = outage + 1 ;
        end
    end
    P_outage(j) = outage/total_syms;      %outage probability
    %here we pass the signal through lowpass filter
    inphase_filtered = filter(b,1,inphase_rx);
    quadphase_filtered = filter(b,1,quadphase_rx);
    
    inphase_value = [];
    quadphase_value = [];
    error(j) = 0;
    for i = 1 : total_syms
        inphase_value(i) = sign(mean(inphase_filtered(tl*(i-1)+1:tl*i)));
        quadphase_value(i) = sign(mean(quadphase_filtered(tl*(i-1)+1:tl*i)));
        %checking error wrt MAP decoder
        if inphase_value(i)~=sign(data3(1,i))
            error(j) = error(j)+1;
        end
        if quadphase_value(i)~=sign(data3(2,i))
            error(j) = error(j) + 1;
        end
    end
    j = j + 1;
end
Pe = error/total_bits;
A = 10^(-3):10^(-3):10^(-2);
N_p = 10^(-13);
SNR_dB = pow2db(A.^2/N_p); %converting amplitude into SNR
figure;
plot(SNR_dB, Pe);axis tight;xlabel('SNR in dB');ylabel('Bit Error rate');grid on;

figure;
plot(SNR_dB,P_outage);xlabel('SNR in dB');ylabel('Outage Probability');grid on;

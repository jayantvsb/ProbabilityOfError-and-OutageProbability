close all;
clear;

fc = 800*10^6;    %carrier frequency
fb = 2*10^6;      %bit frequency is twice of the symbol frequency
Tb = 1/fb;        %bit time period
Th = -135;        %given value of outage threshold
Ts = 1/(4*fc);    %sampling time should be less than the twice of the carrier time period
D = 0.1;          %given distance between tx and rx in km
j = 1;            %iteration number
S = 3;            %no of repeated symbols
%%
total_syms = 1000;      %total no of symbols to be taken
total_bits = 2*total_syms;    %total no of bits
data1 = round(rand(1,total_bits));  %generating random bits of 0 and 1
data2 = 2*data1 - 1;     %converting 0 to -1 and 1 to 1
data3 = reshape(data2,2,length(data1)/2); %here we separate the inphase and quadphase bits
data3_rep = repelem(data3,1,S);         %repeated bits
t = 0:Ts:Tb-Ts; 
tl = length(t);
[tx_bits, Xi, Xq] = modulate(t,fc,data3_rep);    %modulating the transmit bits
L = size(tx_bits,2);
b = fir1(50,0.05,'low');    %here we create low pass filter for demodulation

for A = 10^(-3):10^(-3):10^(-2)
    tx_signal = A*tx_bits;   %amplitude of transmitting signal

    cg_chnl = sqrt(0.5)*randn(1,S*total_syms) + sqrt(0.5)*1i*randn(1,S*total_syms);  %here we generate complex gaussian channel based on rayleigh parameter
    h = repelem(cg_chnl,tl);  %here we create slow fading channel with same channel gain for one bit period
    pathloss_dB = 128.1 + 37.6*log10(D);
    pathloss = 10^-(pathloss_dB/20);   %pathloss

    
    P = 10^((-100-30)/10);   %calculation of power
    sigma = sqrt(P/2);
    noise = sigma*randn(1,L)+1i*sigma*randn(1,L);   %randomly generated noise 
    
    signal = pathloss*h.*tx_signal;    %signal without noise
    signal_rx = signal+noise;          %signal with noise

    conj_h = conj(h);
    Yh = real(conj_h.*signal_rx./abs(h));  

    for i = 1 : S*total_syms
        inphase_rx((i-1)*tl+1:i*tl) = Yh((i-1)*tl+1:i*tl).*Xi;
        quadphase_rx((i-1)*tl+1:i*tl) = Yh((i-1)*tl+1:i*tl).*Xq;
    end

    inphase_filtered = filter(b,1,inphase_rx);
    quadphase_filtered = filter(b,1,quadphase_rx);
    
    for i = 1 : S*total_syms
        inphase_rep(i) = mean(inphase_filtered(tl*(i-1)+1:tl*i));
        quadphase_rep(i) = mean(quadphase_filtered(tl*(i-1)+1:tl*i));
        pow_sym = (1/tl)*sum(real(signal((i-1)*tl+1:i*tl)).^2);
        pow_sym_dBm(i) = pow2db(pow_sym)+30;
    end
    inphase_rep = sign(inphase_rep);
    quadphase_rep = sign(quadphase_rep);
    %considering the majority of the receieved bits
    for i = 1 : total_syms
        inphase_value(i) = mode(inphase_rep(S*(i-1)+1:S*i));
        quadphase_value(i) = mode(quadphase_rep(S*(i-1)+1:S*i));
    end
    error(j) = 0;
    outage = 0;
    for i = 1 : total_syms
        if inphase_value(i)~=sign(data3(1,i))
            error(j) = error(j)+1;
        end
        if quadphase_value(i)~=sign(data3(2,i))
            error(j) = error(j) + 1;
        end
        if(sum(pow_sym_dBm(S*(i-1)+1:S*i) < -135) == S)
            outage = outage + 1;
        end
    end
    P_outage(j) = outage/total_syms;
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
function [tx_signal,Xi,Xq] = modulate(t,fc,data)
no_of_bits = size(data,2);
Xi = cos(2*pi*fc*t);
Xq = -sin(2*pi*fc*t);
inphase_sym = [];
quadphase_sym = [];
for i = 1 : no_of_bits
    inphase_sym = [inphase_sym Xi*data(1,i)];
    quadphase_sym = [quadphase_sym Xq*data(2,i)];
end

tx_signal = inphase_sym + quadphase_sym;
end
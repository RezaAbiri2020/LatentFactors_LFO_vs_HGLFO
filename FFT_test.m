
ch=10;
Fi=5;
%data=ECoGdata_2(ch,:);
data=HG_FingersWholeBand(Fi).PureEnv(:,ch);

L=length(data);
wq=fft(data);
P2 = abs(wq/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
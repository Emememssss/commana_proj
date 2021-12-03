%% Amplitude Modulation
%Read the audio signal [Modulating Signal]
[y, Fs] = audioread('epic_voice.wav')
y=y(2:10000,1)
%sound(y, Fs)
figure(1); plot(y)
shg
xlabel 'Time'
ylabel 'Audio Signal'
title('Message Signal')
Fs=100
%%
Fc = 1135000
dt = 1/Fs
t=(0:1:239600)
t=linspace(1,10000)
carrier_signal = sin(2*pi*t)
figure(2); plot(carrier_signal)
xlabel('Time')
ylabel('Carrier Amplitude')
title('Carrier Signal')
shg
%%
% modulated_signal = (y+0.6).*carrier_signal
%If modulation index is assumed one (1) then Am and Ac is 0.6
Modulated = y.*carrier_signal
figure(3)
plot(t,Modulated)

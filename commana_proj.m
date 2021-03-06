%% Amplitude Modulation
%Read the audio signal [Modulating Signal]
[y, Fs] = audioread('epic_voice.wav')
%sound(y, Fs)
figure(1); plot(y)
shg
xlabel 'Time'
ylabel 'Audio Signal'
title('Message Signal')
%%
Fc = 1135000
t=linspace(1,10000)
n=[0:1/10:30000]
carrier_signal = sin(2*pi*n)
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
modulated_r = Modulated'
plot(t, modulated_r)

%% pre test
n = [0:1/10:2]
y = sin(2*pi*n)
plot(y)

%% test

minLength = min(length(y), length(carrier_signal))
y = y(1:minLength)
carrier_signal = carrier_signal(1:minLength)

%% test 2
%Am =(Vc+Vm*sin(2*pi*fm*t)).*sin(2*pi*fc*t)
Am = (1+y).*carrier_signal
figure(3)
plot(Am)
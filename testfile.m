clear all

%% close all

%% Modulating Signal
m = 1

Vm = 20
fm = 150000
Tm = 1/fm
t = 0:Tm/999:6*Tm
vm = Vm*sin(2*pi*fm*t)
% figure(1);
% subplot(3,1,1);
figure(1);
plot(t,vm);
title('Modulating Signal')

%% Carrier Signal

Vc=Vm/m
fc = 10000000
Tc = 1/fc
vc = Vc*sin(2*pi*fc*t)
% figure(1);
% subplot(3,1,2);
figure(2);
plot(t,vc)
title('Carrier Signal')

%% Amplitude Modulation

Am = (Vc+Vm*sin(2*pi*fm*t)).* sin(2*pi*fc*t);
% subplot(3,1,3);
figure(3);
plot(t,Am);
title('Amplitude Modulation')

%% Project Proper || Amplitude Modulation
test = vc.*vm
figure(10);
plot(test)

%% Modulated Signal
%Read the audio signal [Modulating Signal]
[voice_message_signal, Fs] = audioread('epic_voice.wav')
%sound(y, Fs)
figure(4); plot(voice_message_signal)
shg
xlabel 'Time'
ylabel 'Audio Signal'
title('Message Signal')

%% Carrier Signal

Vc =Vm/m
fc = 10000000
Tc = 1/fc
n = 0:Tm/999:300*Tm
vc = Vc*sin(2*pi*fc*n)
% figure(1);
% subplot(3,1,2);
figure(2);
plot(n,vc)
title('Carrier Signal')

%% Length Adjustment

carrier_signal = vc

minLength = min(length(voice_message_signal), length(vc))
% minLength = minLength/2
voice_message_signal = voice_message_signal(1:minLength);
vc = vc(1:minLength);


%% Voice Modulated Signal
Voice_Modulated_Signal = voice_message_signal .* vc
% Voice_Modulated_Signal = (Vc+voice_message_signal) .* sin(2*pi*fc*t)
figure(5);
plot(minLength,Voice_Modulated_Signal)

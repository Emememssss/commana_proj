clear all
close all

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
voice_message_signal = voice_message_signal(1:minLength)
vc = vc(1:minLength)


%% Voice Modulated Signal
Voice_Modulated_Signal = voice_message_signal .* vc
% % Voice_Modulated_Signal = (Vc+voice_message_signal) .* sin(2*pi*fc*t)
figure(5);
plot(minLength,Voice_Modulated_Signal)

%% AM Receiver

% RF section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fc = [700 750 800]*1e3;     % Carrier frequencies (Hz)
% Ac = [1.00 1.25 1.50];      % Carrier amplitudes
% 
% Fm = [1 2 3]*1e3;           % Modulation frequencies (Hz)
% Dm = [0.25 0.25 0.25];      % Modulation depths
% 
% Fs = 20*max(Fc);            % Sample rate, 20 times the highest RF (Hz)
Ts = 1/Fs;                  % Sample period (s)
% 
% L = 10/min(Fm);             % Duration of signal, 10 times the period of
%                             % the lowest modulation frequency
% 
% t = Ts*(0:ceil(L/Ts)-1);    % Array of sample times (s)
% 
% Sc = diag(Ac)*cos(2*pi*Fc'*t);  % Carrier signals. A three row array with
%                                 % each row representing a single RF
%                                 % carrier.
% 
% Sm = 1 + diag(Dm)*cos(2*pi*Fm'*t);  % Modulating signals. A three row array
%                                     % with each row representing the
%                                     % modulation for a single carrier.
% 
% Stx = sum(Sm.*Sc, 1);   % RF signal. The superposition of three separately
%                         % modulated carriers. This is the type of signal
%                         % that could be expected at the output of the LNA
%                         % (or input to the mixer).

% minLength = min(length(Am), length(vc))
% Am = Am(1:minLength)
% vc = vc(1:minLength)
% vc = vc(1:5995)
% d = Am .* vc
% d1 = conv(d,exp(-t/fm))
% 
% figure(999); plot(t(1:1000),d1(1:1000))

% % zzzz = (1/pi)*(ac+mt)
% zzzz = (1/pi)*(Vc+vm)
% figure(999); plot(zzzz)


zzzz = (1/pi)*(Vc+Voice_Modulated_Signal)
figure(999); plot(zzzz)
sound(zzzz,Fs)

test = lowpass(Voice_Modulated_Signal,2*pi*Fc)

%% Mixer section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% FLO = 300e3;                % Local oscillator frequency (Hz)
% ALO = 1;                    % Local oscillator amplitude
% SLO = ALO*cos(2*pi*FLO*t);  % Local oscillator signal
% 
% Smix = Stx.*SLO;            % Signal at the output of the mixer

%% IF filter section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% There is probably a slicker way to write this in MATLAB, but I am an
% analog guy and don't really have a good handle on the nuances of digital
% filtering. So I just generated a continous time transfer function for a
% Butterworth bandpass filter and then coverted that to its discrete
% equivalent.
% 
% [NUM, DEN] = butter(5, [2*pi*430e3 2*pi*470e3], 's'); % Filter coefficients
%                                                       % for a 10th order
%                                                       % Butterworth
%                                                       % bandpass centered
%                                                       % at 450 MHz
% 
% [NUM, DEN] = butter(5, [2*pi*990000 2*pi*1100000], 's');
% 
% Hd = c2d(tf(NUM, DEN), Ts);     % Discrete equivalent derived from previous
%                                 % continuous time filter coefficients
% 
% Sfilt = filter(Hd.num{1}, Hd.den{1}, Smix);     % Signal at the output of
%                                                 % the IF filter

%% Envelope detector section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Srect = Sfilt; Srect(Srect<0) = 0;  % Half-wave rectified IF signal
% 
% tau = 0.1e-3;                       % Filter time constant (s)
% 
% a = Ts/tau;
% Srect_low = filter(a, [1 a-1], Srect);  % Low pass filtering to recover the
%                                         % the modulating signal

%% Plotting section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                        
                                        
% The plots display numerical data from somewhere in middle of the arrays
% so that the transient responses from the filters have had a chance to
% ring out. Each figure contains three plots: the RF signal, the IF filter
% output, and the demodulated audio signal. The first figure plots a longer
% segment of time so the demodulated audio signal can be distinguished.
% The second figure plots a much shorter segment of time to show the detail
% in the RF signal.
% 
% figure;
% min_index = ceil(length(t)/2);
% max_index = min_index + ceil(2/min(Fm)/Ts);
% 
% subplot(3,1,1);
% plot(t(min_index:max_index), Stx((min_index:max_index)));
% xlim([t(min_index) t(max_index)]); xlabel('Time (s)');
% 
% subplot(3,1,2);
% plot(t(min_index:max_index), Sfilt((min_index:max_index)));
% xlim([t(min_index) t(max_index)]); xlabel('Time (s)');
% 
% subplot(3,1,3);
% plot(t(min_index:max_index), Srect_low((min_index:max_index)));
% xlim([t(min_index) t(max_index)]); xlabel('Time (s)');
% 
% figure;
% min_index = ceil(length(t)/2);
% max_index = min_index + ceil(150/min(Fc)/Ts);
% 
% subplot(3,1,1);
% plot(t(min_index:max_index), Stx((min_index:max_index)));
% xlim([t(min_index) t(max_index)]); xlabel('Time (s)');
% 
% subplot(3,1,2);
% plot(t(min_index:max_index), Sfilt((min_index:max_index)));
% xlim([t(min_index) t(max_index)]); xlabel('Time (s)');
% 
% subplot(3,1,3);
% plot(t(min_index:max_index), Srect_low((min_index:max_index)));
% xlim([t(min_index) t(max_index)]); xlabel('Time (s)');

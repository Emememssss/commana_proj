
[y, Fs] = audioread('epic_voice.wav')
sound(y, Fs)
plot(y)
xlabel 'Time'
ylabel 'Audio Signal'
shg
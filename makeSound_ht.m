%% Create sound: harmonic tone
% Copyright (c) 2019, Sijia Zhao.  All rights reserved.
% Contact: sijia.zhao.10@ucl.ac.uk

clear;

Fs = 44100;
stimDur = 0.5;       % in [sec]
freqStep = 200;
nonComponent = 30;

Ns  = floor(Fs*stimDur);
t  = 1/Fs:1/Fs:stimDur;
nr = 150;
R  = sin(linspace(0,pi/2,nr));
R  = [R,ones(1,Ns-nr*2),fliplr(R)];

y = zeros(size(R));
for a = 1:30
    freq = a*200;
    t0 = floor(Fs/2*rand);
    y = y + sin(2*pi*freq*(t+t0));
end
signal = y.*R;
signal = signal/max(abs(signal));
signal = signal*0.8;

signal = framp(30,signal,Fs); % add 30ms ramp on both sides

filename = ['ht_' num2str(freqStep) 'Hz.wav'];
audiowrite(filename,signal,Fs);

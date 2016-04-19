% HW3 Sample codes 
% Pseudo-Parking Sonar
%
%                                               Edited by Meng-Lin Li, 04/07/2016

% Linear FM/Chirp signal as shown in Fig. 6.15(a) of the textbook
close all
Fs = 100; % Sampling rate in Hz
tau = 10; % Time duration of the linear FM signal in sec
B = 10; % in Hz
beta = B/tau;
A = 1;
t = 0:1/Fs:tau; % time axis
x = A*sin(pi*beta*t.^2); % the linear FM signal
x = hamming(length(x)).'.*x; % hamming windowed linear FM signal
f_axis = 0:1/tau:(length(x)-1)/tau;
f_axis = f_axis-f_axis(ceil(length(f_axis)/2));
len_f=length(f_axis)
FIR = zeros([1 len_f]);
FIR(1:15*tau)=1;
FIR(end-15*tau:end)=1;
figure
plot(t, x)
xlabel('Time (sec)');
ylabel('Amplitude')
title('Hamming windowed linear FM signal')
figure 
x_fft = fft(x);
length(x_fft)
plot(f_axis,abs(x_fft))
load SonarSignalWithWhiteGausNoise
% y: the noisy echo signal, which is recorded at a sampling rate of Fs (in Hz)

% (a)
% remove the noise and increase the SNR by low pass filtering or band pass filtering (depends on the signal frequency band)
% design the proper filter based on the moving average filter
h = [1 3 3 1]/8;% FIR filter
% the output length L of linear convolution will be equal to M+N-1, where M is the length of y and N is the length of h. 
% !!! You have to take care the extra N-1 points in the filter output in order to correctly estimate the echo time.
%fourier analaysis
sample_num = length(y);
dur= sample_num/Fs;
t_axis=0:1/Fs:(length(y)-1)/Fs;
f_axis = 0:1/dur:(length(y)-1)/dur;
f_axis = f_axis-f_axis(ceil(length(f_axis)/2));
FIR = zeros(size(f_axis));
FIR(1:15*dur)=1;
FIR(end-15*dur:end)=1;
figure
plot(t_axis,y);

% y_fft = fft(y);
% figure
% plot(f_axis,FIR)
% y_fft = y_fft.*FIR;
% y_lpf = abs(ifft(y_fft));
% figure
% plot(xcorr(y_lpf,x));
%fourier analaysis
LPF=1;
for i=1:10
    LPF = conv(LPF, h); % 
end
%fourier analaysis
sample_num = length(y);
dur= sample_num/Fs;
f_axis = 0:1/dur:(length(y)-1)/dur;
f_axis = f_axis-f_axis(ceil(length(f_axis)/2));
y_fft = abs(fft(LPF));
%fourier analaysis end
%% result analysis 
y_FIR =conv(y,LPF);
figure
plot(y_FIR)
figure
delay_pdf=xcorr(y_FIR,x);
delay_pdf = delay_pdf((length(delay_pdf)+1)/2:end);
new_taxis = 0:1/Fs:(length(delay_pdf)-1)/Fs;
plot(new_taxis,delay_pdf)
delay_time = t_axis(find(y==max(y)))
save('autocorr.mat','y');
% (b)
% matched filtering
% the output length L of linear convolution will be equal to M+N-1, where M is the length of y and N is the length of x. 
% !!! You have to take care the extra N-1 points in the filter output in order to correctly estimate the echo time.
y_Matched = conv(y, h); 

% (c)
% Check if waveform distortion occurs after filtering, and elaborate the reason

%(d) 
 
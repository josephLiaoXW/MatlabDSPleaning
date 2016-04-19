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

figure
plot(t, x)
xlabel('Time (sec)');
ylabel('Amplitude')
title('Hamming windowed linear FM signal')

figure 
x_fft = fftshift(fft(x));
length(x_fft)
plot(f_axis,abs(x_fft))
xlabel('frquency (sec)');
ylabel('mag')
load SonarSignalWithWhiteGausNoise
match_filter = fliplr(x);
xn_coef = x*x.';
% y: the noisy echo signal, which is recorded at a sampling rate of Fs (in Hz)

% (a)
% remove the noise and increase the SNR by low pass filtering or band pass filtering (depends on the signal frequency band)
% design the proper filter based on the moving average filter
h = [1 3 3 1]/8;% FIR filter
% the output length L of linear convolution will be equal to M+N-1, where M is the length of y and N is the length of h. 
% !!! You have to take care the extra N-1 points in the filter output in order to correctly estimate the echo time.
%fourier analaysis
cutoff_fre = 10;
sample_num = length(y);
dur= sample_num/Fs;
t_axis=0:1/Fs:(length(y)-1)/Fs;
f_axis = 0:1/dur:(length(y)-1)/dur;
f_axis = f_axis-f_axis(ceil(length(f_axis)/2));
FIR = zeros(size(f_axis));
zerof_pos= ceil(length(f_axis)/2);
FIR(zerof_pos-cutoff_fre*dur:zerof_pos+cutoff_fre*dur)=1;
%fourier analysis for y
y_fft = fftshift(fft(y));

figure
plot(f_axis,abs(y_fft));
title('Fourier analysis y')
xlabel('frquency (Hz)');
ylabel('mag')

y_lpf=y_fft.*FIR;

figure
plot(t_axis,abs(y_lpf));
title('Fourier analysis y lpf')
xlabel('frquency (Hz)');
ylabel('mag')

y_lpf = abs(ifft(y_lpf));

figure
plot(t_axis,y_lpf);
xlabel('Time (sec)');
ylabel('Amplitude')
title('y lpf signal')

delay_lpf=xcorr(y_lpf,x);
delay_lpf = delay_lpf((length(delay_lpf)+1)/2:end);

figure
plot(t_axis,delay_lpf);
title('Fourier LPF reslut')
xlabel('Time (sec)');
ylabel('auto-corr')

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
yn_coeff = conv(y.*y,ones(size(x)));
yn_coeff = yn_coeff(length(x):end);
y_FIR =conv(y,LPF);

figure
delay_pdf=xcorr(y_FIR,x);
delay_pdf = delay_pdf((length(delay_pdf)+1)/2:end);
new_taxis = 0:1/Fs:(length(delay_pdf)-1)/Fs;
plot(new_taxis,delay_pdf)

delay_time = t_axis(find(delay_pdf==max(delay_pdf)))

save('autocorr.mat','y');

% delay_pdf=xcorr(y,x);
% delay_pdf = delay_pdf((length(delay_pdf)+1)/2:end);
% figure
% plot(t_axis,delay_pdf);
% title('acorr Filter')
% xlabel('Time (sec)');
% ylabel('auto-corr')

% (b)
% matched filtering
% the output length L of linear convolution will be equal to M+N-1, where M is the length of y and N is the length of x. 
% !!! You have to take care the extra N-1 points in the filter output in order to correctly estimate the echo time.
y_Matched = conv(y, match_filter ); 
y_Matched=y_Matched(length(match_filter):end) ; %%remove Trasient part
y_Matched = y_Matched./((yn_coeff*xn_coef).^(1/2)); %% normalized the mathed signal
figure
plot(y_Matched);
title('Match Filter')
xlabel('Time (sec)');
ylabel('amplitude')
delay_time = t_axis(find(y_Matched==max(y_Matched)))
% (c)
% Check if waveform distortion occurs after filtering, and elaborate the reason

%(d) 
 
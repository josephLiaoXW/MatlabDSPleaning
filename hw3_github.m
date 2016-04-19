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
xf_axis = f_axis-f_axis(ceil(length(f_axis)/2));

figure
plot(t, x)
xlabel('Time (sec)');
ylabel('Amplitude')
title('Hamming windowed linear FM signal')

x_fft = fftshift(fft(x));
figure
subplot(2,1,1)
plot(xf_axis,abs(x_fft)/length(x));
title('Weighted linear FM Fourier analysis')
xlabel('Frequency (Hz)');
ylabel('amplitude')
subplot(2,1,2)
[p_x,w] = phasez(x,length(x));
plot((w-max(w)/2)*100/max(w),fftshift(p_x));
xlabel('angle');
ylabel('amplitude')
load SonarSignalWithWhiteGausNoise
match_filter = fliplr(x);
xn_coef = x*x.';
% y: the noisy echo signal, which is recorded at a sampling rate of Fs (in Hz)

% (a)
% remove the noise and increase the SNR by low pass filtering or band pass filtering (depends on the signal frequency band)
% design the proper filter based on the moving average filter
;% FIR filter
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
subplot(2,1,1)
plot(f_axis,abs(y_fft)/length(y));
title('Fourier analysis y')
xlabel('frquency (Hz)');
ylabel('mag')
subplot(2,1,2)
[p_y,w] = phasez(y,length(y));
plot((w-max(w)/2)*100/max(w),fftshift(p_y));
xlabel('angle');
ylabel('amplitude')


f_axis = 0:1/dur:(length(y)-1)/dur;
f_axis = f_axis-f_axis(ceil(length(f_axis)/2));
%fourier analaysis end
%% result analysis 
yn_coeff = conv(y.*y,ones(size(x)));
yn_coeff = yn_coeff(length(x):end);



% figure
% delay_pdf=xcorr(y_FIR,x);
% delay_pdf = delay_pdf((length(delay_pdf)+1)/2:end);
% new_taxis = 0:1/Fs:(length(delay_pdf)-1)/Fs;
% plot(new_taxis,delay_pdf)
% 
% delay_time = t_axis(find(delay_pdf==max(delay_pdf)))
% 
% save('autocorr.mat','y');

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
delay_time = t_axis(y_Matched==max(y_Matched))

matched_fft = fftshift(fft(y_Matched));

figure
subplot(2,1,1)
plot(f_axis,abs(matched_fft));
title('Match Filter Fourier analysis')
xlabel('Frequency (Hz)');
ylabel('amplitude')
subplot(2,1,2)
[p_matched,w] = phasez(y_Matched,length(y_Matched));
plot((w-max(w)/2)*100/max(w),fftshift(p_matched));
xlabel('angle');
ylabel('amplitude')



ideal_y = zeros(size(y));
ideal_y(delay_time*Fs:delay_time*Fs+length(x)-1) = x;
noise = y_Matched-ideal_y;
noise2 = y-ideal_y;
ps_ideal = ideal_y.*ideal_y;
ideal_match = conv(ideal_y, match_filter ); 
ideal_match=ideal_match(length(match_filter):end) ;
ps_ideal_match = ideal_match.*ideal_match;
psn = noise.*noise;
psn2 = noise2 .*noise2 ;
SNR_org= ps_ideal./psn2;
SNR = ps_ideal_match./psn ;

figure
subplot(2,1,1);
plot(t_axis,ps_ideal)
subplot(2,1,2);
plot(t_axis,psn)

figure
s(1) = subplot(2,1,1);
plot(t_axis,SNR_org)
s(2) = subplot(2,1,2);
plot(t_axis,SNR)
title(s(1),'original SNR')
title(s(2),'filtered SNR')

SNR_org= sum(ps_ideal)./sum(psn2)
SNR = sum(ps_ideal_match)./sum(psn) 

ideal_fft = fftshift(fft(ideal_y));
figure
subplot(2,1,1)
plot(f_axis,abs(ideal_fft));
title('Match Filter Fourier analysis')
xlabel('Frequency (Hz)');
ylabel('amplitude')
subplot(2,1,2)
[p_ideal,w] = phasez(ideal_y,length(ideal_y));
plot((w-max(w)/2)*100/max(w),fftshift(p_ideal));
xlabel('angle');
ylabel('amplitude')
% (c)
% Check if waveform distortion occurs after filtering, and elaborate the reason

%(d) 
 
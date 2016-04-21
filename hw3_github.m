%EE 3660 ¹ù¾ÇÞm u102061222 HW 3 04/20/2016
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
ylabel('magnitude')
subplot(2,1,2)
plot(xf_axis,unwrap(angle(x_fft)));
xlabel('angle');
ylabel('unwrapped angle')

figure
plot(xf_axis,unwrap(angle(x_fft)));

load SonarSignalWithWhiteGausNoise
match_filter = fliplr(x);
xn_coef = x*x.';
% y: the noisy echo signal, which is recorded at a sampling rate of Fs (in Hz)

% (a)
% remove the noise and increase the SNR by low pass filtering or band pass filtering (depends on the signal frequency band)
% design the proper filter based on the moving average filter
% FIR filter
% the output length L of linear convolution will be equal to M+N-1, where M is the length of y and N is the length of h. 
% !!! You have to take care the extra N-1 points in the filter output in order to correctly estimate the echo time.
%fourier analaysis
sample_num = length(y);
dur= sample_num/Fs;
t_axis=0:1/Fs:(length(y)-1)/Fs;
f_axis = 0:1/dur:(length(y)-1)/dur;
f_axis = f_axis-f_axis(ceil(length(f_axis)/2));
%fourier analysis for y
y_fft = fftshift(fft(y));

figure
subplot(2,1,1)
plot(f_axis,abs(y_fft)/length(y));
title('Fourier analysis y')
xlabel('frquency (Hz)');
ylabel('mag')
subplot(2,1,2)
plot(f_axis,unwrap(angle(y_fft)));
xlabel('angle');
ylabel('unwrapped angle')



f_axis = 0:1/dur:(length(y)-1)/dur;
f_axis = f_axis-f_axis(ceil(length(f_axis)/2));
%fourier analaysis end
%% result analysis 
yn_coeff = conv(y.*y,ones(size(x)));
yn_coeff = yn_coeff(length(x):end);



cutoff = 8;
FIR=[zeros(1,(length(y)-cutoff)/2),ones(1,cutoff),zeros(1,(length(y)-cutoff)/2)]/cutoff;


y_FIR = conv(FIR,y);
y_FIR = y_FIR((length(FIR)-cutoff)/2:(length(FIR)-cutoff)/2+length(y)-1);
k=fftshift(fft(FIR));
figure
plot(f_axis,unwrap(angle(k)));
figure
a(1)=subplot(4,1,1);
plot(f_axis,abs(k));
a(4)=subplot(4,1,2);
plot(f_axis,unwrap(angle(k)));
y_FIR = conv(FIR,y);
y_FIR = y_FIR((length(FIR)-cutoff)/2:(length(FIR)-cutoff)/2+length(y)-1);
a(2)=subplot(4,1,3);
plot(t_axis,y_FIR)
k=fftshift(fft(y_FIR));
a(3)=subplot(4,1,4);
plot(f_axis,abs(k));
title(a(1),['Frequency response of FIR (magnitude and phase)'])
title(a(2),['filtered y'])
title(a(3),['Frequency response of filtered y'])
xlabel(a(1),'frquency (Hz)');
ylabel(a(1),'mag')
xlabel(a(4),'frquency (Hz)');
ylabel(a(4),'unwrapped angle')
xlabel(a(2),'time (t)');
ylabel(a(2),'amplitude')
xlabel(a(3),'frquency (Hz)');
ylabel(a(3),'mag')
%Assume noise will not distortion signal too much: the max output is the
%same position as the max amplitude of delay x
%calculate the delay_time
[a,vv] = max(y_FIR); %%postion of max y FIR is the same as max of x
delay_time_FIR = (vv-ceil(length(x)/2))./Fs
%SNR for (a)
ideal_y = zeros(size(y));
ideal_y(delay_time_FIR*Fs:delay_time_FIR*Fs+length(x)-1) = x;
noise = y_FIR-ideal_y;
noise2 = y-ideal_y;
ps_ideal = ideal_y.*ideal_y;
ideal_FIR = conv(ideal_y, FIR ); 
ideal_FIR=ideal_FIR((length(FIR)-cutoff)/2:(length(FIR)-cutoff)/2+length(ideal_y)-1) ;
ps_ideal_FIR = ideal_FIR.*ideal_FIR;
psn = noise.*noise;
psn2 = noise2 .*noise2 ;
SNR_org= ps_ideal./psn2;
SNR_FIR= ps_ideal_FIR./psn ;
%show the result
figure
s(1) = subplot(2,1,1);
plot(t_axis,SNR_org)
s(2) = subplot(2,1,2);
plot(t_axis,SNR_FIR)
SNR_org_all= sum(ps_ideal)./sum(psn2)
SNR_all = sum(ps_ideal_FIR)./sum(psn)
title(s(1),['original SNR¡G',num2str(SNR_org_all)])
title(s(2),['filtered SNR¡G',num2str(SNR_all)])
ylabel(s(1),'s(t) energy / n(t) energy')
ylabel(s(2),'s(t) energy / n(t) energy')
xlabel(s(1),'time (s)')
xlabel(s(2),'time (s)')
% (b)
% matched filtering
% the output length L of linear convolution will be equal to M+N-1, where M is the length of y and N is the length of x. 
% !!! You have to take care the extra N-1 points in the filter output in order to correctly estimate the echo time.


y_Matched = conv(y, match_filter ); 
y_Matched=y_Matched(length(match_filter):end) ; %%remove Trasient part
y_Matched = y_Matched./((yn_coeff*xn_coef).^(1/2)); %% normalized the mathed signal
figure
plot(t_axis,y_Matched);
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
ylabel('mag')
subplot(2,1,2)
[p_matched,w] = phasez(y_Matched,length(y_Matched));
plot(f_axis,unwrap(angle(matched_fft)));
xlabel('angle');
ylabel('unwrapped angle')



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
SNR_org_all= sum(ps_ideal)./sum(psn2)
SNR_all = sum(ps_ideal_match)./sum(psn)
title(s(1),['original SNR¡G',num2str(SNR_org_all)])
title(s(2),['filtered SNR¡G',num2str(SNR_all)])

 

ideal_fft = fftshift(fft(ideal_y));
figure
subplot(2,1,1)
plot(f_axis,abs(ideal_fft));
title('Ideal X delay without noise')
xlabel('Frequency (Hz)');
ylabel('Magnitude')
subplot(2,1,2)
plot(f_axis,unwrap(angle(ideal_fft)));
xlabel('angle');
ylabel('unwrapped angle')
%(c)
%elaborate the distortion of x through matched fiter
figure
subplot(3,1,1)
plot(t_axis,ideal_match);
title('distortion of match filter')
xlabel('time (s)');
ylabel('amplitude')
subplot(3,1,2)
plot(f_axis,abs(fftshift(fft(ideal_match))));
xlabel('Frequency (Hz)');
ylabel('magnitude')
subplot(3,1,3)
[p_ideal,w] = phasez(ideal_y,length(ideal_match));
plot(f_axis,unwrap(angle(fftshift(fft(ideal_match)))));
xlabel('angle');
ylabel('unwrapped angle')
%(d)
k=10;
y_intmatch = interp(y_Matched,k);
t_axis=0:1/Fs/k:(length(y_intmatch)-1)/Fs/k;
delay_time = t_axis(y_intmatch==max(y_intmatch))
k=10;
y_interpo = interp(y,k);
x_interpo = interp(x,k);
Fs=Fs*k;
t_axis=0:1/Fs:(length(y_interpo)-1)/Fs;

match_filter = fliplr(x_interpo);
xn_coef = x_interpo*x_interpo.';
yn_coeff = conv(y_interpo.*y_interpo,ones(size(x_interpo)));
yn_coeff = yn_coeff(length(x_interpo):end);

y_iMatched = conv(y_interpo, match_filter ); 
y_iMatched=y_iMatched(length(match_filter):end) ; %%remove Trasient part
y_iMatched = y_iMatched./((yn_coeff*xn_coef).^(1/2)); %% normalized the mathed signal
figure
plot(t_axis,y_iMatched);
title('Match Filter')
xlabel('Time (sec)');
ylabel('amplitude')
[a,vv] = max(y_iMatched);
delay_time = t_axis(vv)
normalized_time = (vv-1)/k
%(d)(e)
%section for estimating tolerance
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
y_d = downsample(y,100/Fs);
match_filter = fliplr(x);
xn_coef = x*x.';
yn_coeff = conv(y_d.*y_d,ones(size(x)));
yn_coeff = yn_coeff(length(x):end);


y_iMatched = conv(y_d, match_filter ); 
y_iMatched=y_iMatched(length(match_filter):end) ; %%remove Trasient part
y_iMatched = y_iMatched./((yn_coeff*xn_coef).^(1/2)); %% normalized the mathed signal
t_axis=0:1/Fs:(length(y_d)-1)/Fs;
k=5;
y_intmatch = interp(y_iMatched,k,2,12/50);
t_axis=0:1/Fs/k:(length(y_intmatch)-1)/Fs/k;
[a,vv] = max(y_intmatch);
delay_time = t_axis(vv)
normalized_time = (vv)/(Fs*k/100)

ideal_y = zeros(size(y_d));
ideal_y(round(delay_time*Fs):round(delay_time*Fs)+length(x)-1) = x;
y_iMatched = [y_iMatched,zeros([1 length(ideal_y)-length(y_iMatched)])];
noise = y_iMatched-ideal_y;
ps_ideal = ideal_y.*ideal_y;
ideal_match = conv(ideal_y, match_filter ); 
ideal_match=ideal_match(length(match_filter):end) ;
ps_ideal_match = ideal_match.*ideal_match;
psn = noise.*noise;
SNR = ps_ideal_match./psn;
SNR_all = sum(ps_ideal_match)./sum(psn)
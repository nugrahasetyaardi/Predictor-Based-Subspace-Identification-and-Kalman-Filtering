function [y_filtered] = fft_and_filter(y_0,dt,fc,order,name)
warning off
% fft_and_filter gives the single-sided spectrum and the filtered signal
%
% The inputs are:
% > the discrete signal
% > sampling time
% > cutoff frequency
% > order of the filter
% > name of the signal
%
% The output is the filtered signal

%% FFT computation

L = length(y_0);
y_fft = fft(y_0);
Fs = 1/dt; % sampling frequency
f = Fs*(0:(L/2))/L;% generate frequency vector for FFT plot

% Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
P2 = abs(y_fft/L);
P1= P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure
subplot(2,1,1);
plot(f,P1);
title('Single-Sided Spectrum')
xlabel('f (hz)')
ylabel('|P1(f)|')
grid on;grid minor;

%% Butter filter

%zero-phase digital filtering by processing the input data, y_0
%in both the forward and reverse directions

[b,a] = butter(order,fc/(Fs/2));
y_filtered = filtfilt(b,a,y_0);

subplot(2,1,2);
plot(y_0)
hold on
plot(y_filtered,'black');
legend('Measured','Filtered');
title(name);xlabel('Samples'); ylabel(name);grid on;


end


%% Open data 
open("Nikolai hard force export.csv") % open .csv file with importer, select "column vector" for data type
%% Plot stuff
EMG1 = rmmissing(Ch1)
time = 11.3595 %Insert final time in seconds here   //Insert final time from the excel sheet
d = length(EMG1)
ElapsedTime = linspace(0,time,d)
figure(1)
plot(ElapsedTime,EMG1)
title('EMG data')
ylabel('EMG(mV)')
xlabel('Time(s)')

%% Pass Filters - Step 2

% [B,A] = butter(N,Wn,'low') designs a highpass filter where all
% frequencies below the value are retained and the ones above this
% value are filtered
Fs = 1/(ElapsedTime(2)-ElapsedTime(1));
Fcut = 500; % Filter cut-off frequency
N = 2; % Filter order (for filtfilt it becomes a 4th order)
Wn = (Fcut)/(Fs/2); % Normalized cut off frequency
[b,a] = butter(N,Wn,'low');           % IIR filter design

% Apply the filter
a = filtfilt(b,a,EMG1); % zero-phase filtering

figure(2)
plot(ElapsedTime,a)
title('Low Pass Filter')
ylabel('EMG(mV)')
xlabel('Time(s)')
%% Pass Filters - Step 3

% [B,A] = butter(N,Wn,'low') designs a highpass filter where all
% frequencies below the value are retained and the ones above this
% value are filtered

Fcut2 = 10; % Filter cut-off frequency
N = 2; % Filter order (for filtfilt it becomes a 4th order)
Wn = (Fcut2)/(Fs/2); % Normalized cut off frequency
[d,c] = butter(N,Wn,'high');           % IIR filter design

% Apply the filter
z = filtfilt(d,c,a); % zero-phase filtering

figure(3)
plot(ElapsedTime,z)
title('High Pass Filter')
ylabel('EMG(mV)')
xlabel('Time(s)')
%% Full Wave Rectify - Step 4

absZ = abs(z) 
figure(4) 
plot(ElapsedTime,absZ) 
title('Full Wave Rectify')
ylabel('EMG(mV)')
xlabel('Time(s)')
%% Linear Envelope - Low pass, Butter - Step 5

% [B,A] = butter(N,Wn,'low') designs a highpass filter where all
% frequencies below the value are retained and the ones above this
% value are filtered

Fcut3 = 5; % Filter cut-off frequency
N = 2; % Filter order (for filtfilt it becomes a 4th order)
Wn = (Fcut3)/(Fs/2); % Normalized cut off frequency
[f,e] = butter(N,Wn,'low');           % IIR filter design

% Apply the filter
x = filtfilt(f,e,absZ); % zero-phase filtering

figure(4)
plot(ElapsedTime,x)
title('Linear Envelope')
ylabel('EMG(mV)')
xlabel('Time(s)')
%% FFT - No filters

fftNF = fft(EMG1)
figure(5)
plot(ElapsedTime,fftNF)
title('EMG Data - Hard grip, FFT without filters')
ylabel('EMG(mV)')
xlabel('Time(s)')
%% %% FFT - After filters

fftF = fft(z)
figure(5)
plot(ElapsedTime,fftF)
title('EMG Data - Hard grip, FFT after filters')
ylabel('EMG(mV)')
xlabel('Time(s)')


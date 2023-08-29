clc
clear all
close all

%!lskdb -p 1ms -s 14/06/2023-13:55:19 -e 14/06/2023-14:10:21 -d Kite_pos Kite_control 3axis_FS -m Dynamic


dataRaw = readtable("14062023_135519.csv");

data = dataRaw(:,{'x_time', ...
    'KITE_elevation_wr_Dynamic', ...
    'KITE_azimuth_wr_Dynamic', ...
    'KITE_lacet_wr_Dynamic',...
    'CTRL_dyn_command_Dynamic'});

data = renamevars(data,["x_time", ...
    "KITE_elevation_wr_Dynamic", ...
    "KITE_azimuth_wr_Dynamic", ...
    "KITE_lacet_wr_Dynamic", ...
    "CTRL_dyn_command_Dynamic"],...
    ["time","theta","phi","psi","delta"]);

data(any(ismissing(data),2), :) = [];

rows = height(data);

newTime = rowfun(@convertStringToSeconds, data,"InputVariables","time","OutputVariableNames","time");
data.time = [];
data = [newTime  data];
data.time = data.time - data.time(1);

writetable(data, "data.csv");

figure
subplot(211)
hold on
grid on
plot(data.time, data.theta,"DisplayName","Élévation")
plot(data.time, data.phi,"DisplayName","Azimuth")
plot(data.time, data.psi,"DisplayName","Lacet")
xlabel('Time (s)');
ylabel('Angle (°)')
legend
title("États kite")
subplot(212)
plot(data.time, data.delta)
grid on
xlabel('Time (s)');
ylabel('Différentiel (mm)')
title("Signal de commande")

%% FFT

Fs = 5; % Sampling frequency (Hz)
t = 0:1/Fs:data.time(end);
nfft = 1024;

u = timeseries(data.delta, data.time);
u = resample(u, t);
y = timeseries(data.psi, data.time);
y = resample(y, t);

U = fft(u.Data,nfft);
U = U(1:nfft/2);
mu = abs(U);

Y = fft(y.Data,nfft);
Y = Y(1:nfft/2);
my = abs(Y);
f = (0:nfft/2-1)*Fs/nfft;

figure
subplot(2,4,[1 2 5 6])
hold on
grid on
plot(u)
plot(y)
legend("Commande","Lacet")
xlabel("Time (s)")

subplot(2,4,[3 4 7 8])
hold on
plot(f,mu);
plot(f,my);
xlabel("Frequency (Hz)")
ylabel("Power")
grid on
legend("Commande","Lacet")

%% Least square identification



%% Converts time from influxDB to time in seconds.
function timeInSeconds = convertStringToSeconds(dateTimeString)
    dateTimeString = split(strrep(dateTimeString, 'Z', ''),"T");
    timeCells = split(dateTimeString{2},":");
    timeInSeconds = str2num(timeCells{1})*3600 + str2num(timeCells{2})*60 + str2num(timeCells{3});
end

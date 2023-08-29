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

save("data_pilote_auto","data");

% Converts time from influxDB to time in seconds.
function timeInSeconds = convertStringToSeconds(dateTimeString)
    dateTimeString = split(strrep(dateTimeString, 'Z', ''),"T");
    timeCells = split(dateTimeString{2},":");
    timeInSeconds = str2num(timeCells{1})*3600 + str2num(timeCells{2})*60 + str2num(timeCells{3});
end

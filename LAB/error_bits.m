clc
clear all
close all

% dateTimeArray = datetime(timestampsArray, 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd HH:mm:ss');

dataRaw = readtable("1.csv");
dataRaw = dataRaw(:,3:end);
dataRaw(any(ismissing(dataRaw),2), :) = [];

newTime = rowfun(@convertStringToSeconds, dataRaw,"InputVariables","x_time","OutputVariableNames","x_time");
dataRaw.x_time = [];
dataRaw = [newTime  dataRaw];
dataRaw.time = dataRaw.x_time - dataRaw.x_time(1);

figure
plot(dataRaw.time, dataRaw.GEN_error_bits)
grid on

% Converts time from influxDB to time in seconds.
function timeInSeconds = convertStringToSeconds(dateTimeString)
    dateTimeString = split(strrep(dateTimeString, 'Z', ''),"T");
    timeCells = split(dateTimeString{2},":");
    timeInSeconds = str2num(timeCells{1})*3600 + str2num(timeCells{2})*60 + str2num(timeCells{3});
end

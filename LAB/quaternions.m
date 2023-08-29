clc
clear all
close all

port = serialportlist;

device = serialport(port,115200);
data = readbinblock(device);

% char(output)
%% Clean the workspace

clc
clear all
close all

%% Parameters

debug = true;

max_r = 5;
fitting_window = 201; % has to bee odd
polynome_degree = 20;
filter_size = 21; % has to bee odd

%% Load data

% load('data/aggressive_0.01.mat');
% data = data(:,2:end);
% 
% t = data(1,:)';
% dt = t(2) - t(1);
% 
% attitude = quat2eul(data(8:11,:)');
% response = data(5:7,:)';
% state = [response, attitude, data(12:17,:)'];
% 
% command = data(18,:);
% response = state(:,3);

load('data/motor.mat');
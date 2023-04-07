%% Clean the workspace

clc
clear all
close all

%% Parameters

system = 1; % 1 = motor, 2 = monocopter, 3 = quadcopter

%% Load data

if system == 1
    table = readtable("data/motor.txt");
    [t, command, response] = parse_table(table);
end
if system == 2
    table = readtable("data/monocopter.txt");
    [t, command, response] = parse_table(table);
end
if system == 3
    table = readtable("data/quadcopter.txt");
    [t, command, response] = parse_table(table);
end


%% Identify relative degree

relative_degree = nan(size(command, 2), size(response, 2));
confidence_level = nan(size(command, 2), size(response, 2));
for i = 1:size(command, 2)
    for j = 1:size(response, 2)
        [relative_degree(i,j), confidence_level(i,j)] = ...
            estimate_relative_degree(t, command(:,i), response(:,j));
    end
end

disp('Relative degree:');
disp(relative_degree);
disp('Confidence level:');
disp(confidence_level);

%% Parse table %%

function [t, command, response] = parse_table(table)
    t = table.time;
    command = [];
    c = 1;
    while any(ismember(table.Properties.VariableNames, "command_" + c))
        command = [command, table2array(table(:, "command_" + c))];
        c = c + 1;
    end
    response = [];
    c = 1;
    while any(ismember(table.Properties.VariableNames, "response_" + c))
        response = [response, table2array(table(:, "response_" + c))];
        c = c + 1;
    end
end

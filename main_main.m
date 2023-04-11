%% Clean the workspace

clc
clear all
close all

%% Parameters

system = 3; % 1 = motor, 2 = monocopter, 3 = quadcopter

%% Load data

if system == 1
    table = readtable("data/motor.txt");
    [t, command, response] = parse_table(table);
end

if system == 2
    load('data/monocopter_aggressive_0.01.mat');
    data = data(:,2:end);
    
    t = data(1,:)';
    
    attitude = quat2eul(data(8:11,:)');
    response = data(5:7,:)';
    state = [response, attitude, data(12:17,:)'];
    
    command = data(18,:)';
    response = state(:,1:12);

    % 5  2  0  1  7  0  0  1  1  3  5  5
    response(:,1) = [response(4:end,8);response(end-2:end,8)];
    response(:,2) = response(:,2);
    response(:,3) = [response(1:2,10);response(1:end-2,10)];
    table = array2table([t,command,response(:,1:3)], ...
        'VariableNames', ["time","command_1","response_1","response_2","response_3"]);
    writetable(table, 'data/monocopter.txt');
    writetable(table, 'data/monocopter.csv');

    clear all
    system = 2;

    table = readtable("data/monocopter.txt");
    [t, command, response] = parse_table(table);
end

if system == 3
    table = readtable('data/quadcopter_2.txt');

    attitude = quat2eul([table.quaternion_w table.quaternion_z table.quaternion_y table.quaternion_x]);
    
    t = table.time;
    command(:,1) = table.command_vertical;
    command(:,2) = table.command_roll;
    command(:,3) = table.command_pitch;
    command(:,4) = table.command_yaw;   
    response(:,1) = table.translation_x/1000;
    response(:,2) = table.translation_y/1000;
    response(:,3) = table.translation_z/1000;
    response(:,4) = attitude(:,3);

    ind = 9160:20:numel(t);
    t = t(ind);
    command = command(ind,:);
    response = response(ind,:);

    t = t - t(1);
    response(:,1:3) = response(:,1:3) - response(1,1:3);

    x = response(:,1);
    y = response(:,2);
    response(:,1) = cos(response(:,4)).*x - sin(response(:,4)).*y;
    response(:,2) = sin(response(:,4)).*x + cos(response(:,4)).*y;

    response(5000:end,1) = response(5000,1) + (1:numel(response(5000:end,1)))/1000000;
    response(5000:end,2) = response(5000,2) + (1:numel(response(5000:end,2)))/1000000;

    response(:,4) = response(:,4) - response(1,4);

%         [~, ia, ~] = unique(response, 'rows', 'stable');
%     t = t(ia);
%     command = command(ia,:);
%     response = response(ia,:);

%     command(:,1) = [command(3:end,1);command(end-1:end,1)];
% %     command(:,1) = [command(9:end,1);command(end-7:end,1)];
%     command(:,1) = [command(1:5,1);command(1:end-5,1)];
%     command(:,3) = [command(6:end,3);command(end-4:end,3)];
% %     command(:,3) = [command(1:9,3);command(1:end-9,3)];
% 
%     response(:,1) = [response(2:end,1);response(end-0:end,1)];
    response(:,1) = [response(45:end,1);response(end-43:end,1)];
    response(:,2) = [response(42:end,2);response(end-40:end,2)];
    response(:,3) = [response(4:end,3);response(end-2:end,3)];
    response(:,4) = [response(9:end,4);response(end-7:end,4)];
% %     response(:,4) = [response(1:9,4);response(1:end-9,4)];

    figure
    hold on
    plot(t, command(:,2))
    plot(t, command(:,3))
    plot(t, command(:,4))
    plot(t, response(:,1))
    plot(t, response(:,2))
    xlim([0 90])
    legend('roll', 'pitch', 'yaw', 'x', 'y')

    
%     figure
%     hold on
%     plot(t, command(:,4))
%     plot(t, response(:,4)*100)
%     xlim([100 140])

    table = array2table([t,command,response], ...
        'VariableNames', ["time","command_1","command_2","command_3","command_4","response_1","response_2","response_3","response_4"]);
    writetable(table, 'data/quadcopter.txt');
    writetable(table, 'data/quadcopter.csv');

    clear all
    system = 3;

    table = readtable("data/quadcopter.txt");
    [t, command, response] = parse_table(table);
end


%% Identify relative degree

relative_degree = nan(size(command, 2), size(response, 2));
confidence_level = nan(size(command, 2), size(response, 2));
for i = 1:size(command, 2)
    for j = 1:1%size(response, 2)
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
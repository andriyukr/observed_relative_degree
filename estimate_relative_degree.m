%% Similarity between vectors %%

function [degree, confidence] = estimate_relative_degree(t, command, response)

    %% Parameters
    
    debug = true;
    show_figure = true;
    
    max_r = 9;
    fitting_window = 201; % has to bee odd
    polynome_degree = 20;
    filter_size = 21; % has to bee odd
       
    %% Process data

    dt = t(2) - t(1);
        
    d_command = differentate(command, dt);
    [d_command_peaks.v, d_command_peaks.t] = get_peaks(d_command, t);
    
    %% Interpolate data
    
    for i=1:numel(d_command_peaks.t.min)
        k_min(i) = find(t == d_command_peaks.t.min(i));
        [response_peaks1_v{i}, response_peaks1_t{i}, response_fitted1{i}] = ...
        get_fitted_peaks(k_min(i), t, response, max_r, fitting_window, polynome_degree);   
    end
    for i=1:numel(d_command_peaks.t.max)
        k_max(i) = find(t == d_command_peaks.t.max(i));
        [response_peaks2_v{i}, response_peaks2_t{i}, response_fitted2{i}] = ...
        get_fitted_peaks(k_max(i), t, response, max_r, fitting_window, polynome_degree);   
    end
    
    % for w = 1:fix(numel(t)/100) - 1
    %     window = w*100 - (fitting_window - 1)/2:w*100 + (fitting_window - 1)/2
    %     t_local = t(window) - t(w*100);
    %     command_local = command(window);
    %     responce_local = response(window);
    %     
    %     for i = 1:10
    %         if i == 1
    %             p_response{w,1} = polyfit(t_local, responce_local, 20);
    %         else
    %             p_response{w,i} = polyder(p_response{w,i - 1});
    %         end
    %         response_fitted{w,i} = polyval(p_response{w,i}, t_local(50:150));
    %         [response_peaks_v{w,i}, response_peaks_t{w,i}] = get_peaks(response_fitted{w,i}, t_local(50:150));
    %     end
    % end

    %% Show responce

    if show_figure
        colors = lines(7);

        d_response(:,1) = differentate(response, dt);
        for i = 2:10
            d_response(:,i) = differentate(d_response(:,i - 1), dt);
        end

        figure;
        hold on;
        grid on;
        h_c = plot(t, command/max(abs(command)), 'Color', colors(1,:), 'LineStyle', ':', 'linewidth', 2);
        h_c_d = plot(t, d_command/max(abs(d_command)), 'Color', colors(2,:), 'LineStyle', ':', 'linewidth', 1);
        peaks_t = d_command_peaks.t.min;
        peaks_v = d_command_peaks.v.min/max(abs(d_command));
            scatter(peaks_t, peaks_v, 'MarkerEdgeColor', colors(2,:), 'MarkerFaceColor', [0, 0, 1], 'linewidth', 1);
            peaks_t = d_command_peaks.t.max;
            peaks_v = d_command_peaks.v.max/max(abs(d_command));
            scatter(peaks_t, peaks_v, 'MarkerEdgeColor', colors(2,:), 'MarkerFaceColor', [1, 0, 0], 'linewidth', 1);
        % plot(t(time_window), command_fitted{1}(time_window), 'LineStyle', ':', 'linewidth', 1);
        % plot(t(time_window), command_fitted{2}(time_window), 'LineStyle', ':', 'linewidth', 1);
        h_r = plot(t, response/max(abs(response)), 'Color', colors(3,:), 'linewidth', 2);
        % plot(t, state(:,9), 'linewidth', 2);
        % plot(t, state(:,12), 'linewidth', 2);
        h_r_d1 = plot(t, d_response(:,1)/max(abs(d_response(:,1))), 'Color', colors(4,:), 'linewidth', 1);
        h_r_d2 = plot(t, d_response(:,2)/max(abs(d_response(:,2))), 'Color', colors(5,:), 'linewidth', 1);
        h_r_d3 = plot(t, d_response(:,3)/max(abs(d_response(:,3))), 'Color', colors(6,:), 'linewidth', 1);
        % plot(t, d_response(:,4)/max(abs(d_response(:,4))), 'linewidth', 1);
        for w = 1:numel(k_min) - 1
            window = k_min(w) - (fitting_window - 1)/4:k_min(w) + (fitting_window - 1)/4;
        
            h_f1 = plot(t(window), response_fitted1{w}{1}/max(abs(response)), 'Color', colors(2,:), 'LineStyle', '--', 'linewidth', 2);
            h_f2 = plot(t(window), response_fitted1{w}{2}/max(abs(d_response(:,1))), 'Color', colors(4,:), 'LineStyle', '--', 'linewidth', 2);
            h_f3 = plot(t(window), response_fitted1{w}{3}/max(abs(d_response(:,2))), 'Color', colors(5,:), 'LineStyle', '--', 'linewidth', 2);
            h_f4 = plot(t(window), response_fitted1{w}{4}/max(abs(d_response(:,3))), 'Color', colors(6,:), 'LineStyle', '--', 'linewidth', 2);
            h_f5 = plot(t(window), response_fitted1{w}{5}/max(abs(d_response(:,4))), 'Color', colors(7,:), 'LineStyle', '--', 'linewidth', 2);
        
            for i = 2:5
                peaks_t = response_peaks1_t{w}{i}.min + t(k_min(w));
                peaks_v = response_peaks1_v{w}{i}.min/max(abs(d_response(:,i - 1)));
                scatter(peaks_t, peaks_v, 'MarkerEdgeColor', colors(i + 2,:), 'MarkerFaceColor', [0, 0, 1], 'linewidth', 2);
                peaks_t = response_peaks1_t{w}{i}.max + t(k_min(w));
                peaks_v = response_peaks1_v{w}{i}.max/max(abs(d_response(:,i - 1)));
                scatter(peaks_t, peaks_v, 'MarkerEdgeColor', colors(i + 2,:), 'MarkerFaceColor', [1, 0, 0], 'linewidth', 2);
            end
        end
        for w = 2:numel(k_max)
            window = k_max(w) - (fitting_window - 1)/4:k_max(w) + (fitting_window - 1)/4;
        
            h_f1 = plot(t(window), response_fitted2{w}{1}/max(abs(response)), 'Color', colors(2,:), 'LineStyle', '--', 'linewidth', 2);
            h_f2 = plot(t(window), response_fitted2{w}{2}/max(abs(d_response(:,1))), 'Color', colors(4,:), 'LineStyle', '--', 'linewidth', 2);
            h_f3 = plot(t(window), response_fitted2{w}{3}/max(abs(d_response(:,2))), 'Color', colors(5,:), 'LineStyle', '--', 'linewidth', 2);
            h_f4 = plot(t(window), response_fitted2{w}{4}/max(abs(d_response(:,3))), 'Color', colors(6,:), 'LineStyle', '--', 'linewidth', 2);
            h_f5 = plot(t(window), response_fitted2{w}{5}/max(abs(d_response(:,4))), 'Color', colors(7,:), 'LineStyle', '--', 'linewidth', 2);
        
            for i = 2:5
                peaks_t = response_peaks2_t{w}{i}.min + t(k_max(w));
                peaks_v = response_peaks2_v{w}{i}.min/max(abs(d_response(:,i - 1)));
                scatter(peaks_t, peaks_v, 'MarkerEdgeColor', colors(i + 2,:), 'MarkerFaceColor', [0, 0, 1], 'linewidth', 2);
                peaks_t = response_peaks2_t{w}{i}.max + t(k_max(w));
                peaks_v = response_peaks2_v{w}{i}.max/max(abs(d_response(:,i - 1)));
                scatter(peaks_t, peaks_v, 'MarkerEdgeColor', colors(i + 2,:), 'MarkerFaceColor', [1, 0, 0], 'linewidth', 2);
            end
        end
        % xlim([time_window(1)*dt time_window(end)*dt]);
        % ylim([-100, 100]);
%         xlim([0.19 0.29])
%         ylim([-0.5 0.5])
        
        legend( ...
            [h_c, h_c_d, ...
            h_r, ...
            h_r_d1, h_r_d2, h_r_d3, ...
            h_f1, h_f2, h_f3, h_f4, h_f5], ...
            'input', 'd\_input', ...%'input\_fitted', 'd\_input\_fitted', ...
            'response', ...%'state1', 'state2', ...
            'd\_response\_1', 'd\_response\_2', 'd\_response\_3', ...
             'response\_fitted\_0', 'response\_fitted\_1', 'response\_fitted\_2', 'response\_fitted\_3', 'response\_fitted\_4', ...
            'location', 'southeast');
    end
    
    %% Estimate relative degree
    
    for w = 1:numel(k_min)
        for i = 2:max_r
            r_p1_t{i - 1} = response_peaks1_t{w}{i};
        end
        [error_same(w,:), error_opposite(w,:)] = find_closest_excitation(r_p1_t, max_r);
    end
    for w = 1:numel(k_max)
        for i = 2:max_r
            r_p2_t{i - 1} = response_peaks2_t{w}{i};
        end
        [error_opposite(numel(k_min) + w,:), error_same(numel(k_min) + w,:)] = find_closest_excitation(r_p2_t, max_r);
    end
    
    if debug
        error_same
        error_opposite
        mean(error_same)
        mean(error_opposite)
        sum(isnan(error_same))
        sum(isnan(error_opposite))
    end 
    
    if min(mean(error_same,'omitnan')) < min(mean(error_opposite,'omitnan'))
        [error, degree] = min(mean(error_same,'omitnan')); 
    else 
        [error, degree] = min(mean(error_opposite,'omitnan'));
    end
    
    degree = degree - 1;
    
    confidence = 1/(1 + error);

%     cond = degree >= max_r/2 || degree == 0 || confidence <= 0.94;
%     degree(cond) = inf;
%     confidence(cond) = 0;
end

%% Derivative

% https://en.wikipedia.org/wiki/Finite_difference_coefficient
function d_vector = differentate(vector, dt)
    last = numel(vector);
    d_vector(1) = 0;
    d_vector(1:last - 1) = diff(vector)/dt;
    d_vector(last) = 0;

%     d_vector(1) = (vector(2) - vector(1))/dt;
%     d_vector(2:last - 1) = (vector(3:last) - vector(1:last - 2))/(2*dt);
%     d_vector(last) = (vector(last) - vector(last - 1))/dt;
end

%% Get peaks of the derived polynomes

function [polynome_peaks_v, polynome_peaks_t, polynome_fitted] = get_fitted_peaks(k, time, polynome, max_r, fitting_window, polynome_degree)
    window = k - (fitting_window - 1)/2:k + (fitting_window - 1)/2;
    window = window(window > 0 & window < numel(time));
    time_local = time(window) - time(k);
    polynome_local = polynome(window);

    for i = 1:max_r
        if i == 1
            p{1} = polyfit(time_local, polynome_local, polynome_degree);
        else
            p{i} = polyder(p{i - 1});
        end
        time_local_neighbourhood = time_local((fitting_window - 1)/4 + 1:end - (fitting_window - 1)/4);
        polynome_fitted{i} = polyval(p{i}, time_local_neighbourhood);
        [polynome_peaks_v{i}, polynome_peaks_t{i}] = get_peaks(polynome_fitted{i}, time_local_neighbourhood);
    end
end

%% Filters %%

function vector_filtered = median_filter(vector, window)
    vector_filtered = zeros(size(vector));
    range = (window - 1)/2;
    for i = range + 1:numel(vector) - range
        vector_filtered(i) = median(vector(i - range:i + range));
    end
    vector_filtered(1:range) = vector_filtered(range + 1);
    vector_filtered(end - range:end) = vector_filtered(end - range);
end

function vector_filtered = mean_filter(vector, window)
    vector_filtered = zeros(size(vector));
    range = (window - 1)/2;
    for i = range + 1:numel(vector) - range
        vector_filtered(i) = mean(vector(i - range:i + range));
    end
    vector_filtered(1:range) = vector_filtered(range + 1);
    vector_filtered(end - range:end) = vector_filtered(end - range);
end

%% Find peaks %%

function peaks_positions = find_peaks(vector, n)
    vector_sorted = sort(abs(vector),'descend');
    vector_sorted(1:n);
    i = 1;
    while i < n + 1
        vector_sorted(i);
        n_peaks = numel(find(abs(vector) == vector_sorted(i)));
        n_peaks_max = min([n + 1 - i, n_peaks]);
        peaks_positions(i:i + n_peaks_max - 1) = find(abs(vector) == vector_sorted(i), n_peaks_max);
        i = i + n_peaks_max;
    end
end

function [peaks_p, peaks_t] = get_peaks(p, t)
    peaks_p.min = p(islocalmin(p))';
    peaks_p.max = p(islocalmax(p))';
    peaks_t.min = t(islocalmin(p))';
    peaks_t.max = t(islocalmax(p))';
end

%% Find excitation

function degree = find_first_excitation(peaks_response, peaks_command)
    for i = 1:9
        min_max = [peaks_response{i}.min() peaks_response{i}.max];
        if numel(min_max(min_max >= [peaks_command.min peaks_command.max])) > 0
            excitation(i) = min(min_max(min_max >= [peaks_command.min peaks_command.max]));
        else
            excitation(i) = nan;
        end
    end
    [~, degree] = min(excitation);
end

function [error_min, error_max] = find_closest_excitation(peaks_response, max_r)
    for i = 1:max_r - 1
        if numel(peaks_response{i}.min) > 0 % min
            error_min(i) = min(abs(peaks_response{i}.min));
        else
            error_min(i) = nan;
        end
        if numel(peaks_response{i}.max) > 0 % max
            error_max(i) = min(abs(peaks_response{i}.max));
        % none
        else
            error_max(i) = nan;
        end
    end
end

%% Similarity between vectors %%

function s = similarity1(v1, v2)
    s = mean(ismember(v1, v2));
end

function s = similarity(v1, v2)
    if numel(v1) > 0 && numel(v2) > 0
        for i = 1:numel(v1)
            d1(i) = min(abs(v2 - v1(i)));
        end
        for i = 1:numel(v2)
            d2(i) = min(abs(v1 - v2(i)));
        end
        s = mean(d1) + mean(d2);
    else
        s = nan;
    end
end

%% Remove duplicates %%


function [t, command, response] = remove_duplicates(t, command, response)
    [~, ia, ~] = unique([command response], 'rows', 'stable');
    t = t(ia);
    command = command(ia);
    response = response(ia);
end
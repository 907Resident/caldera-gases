function [fig03,fig04,fig05] = basic_gas_viz(ChamON_data, nchams,       ...
                                             site_tag, ddmmmyyyy,       ...
                                             working_dir)
%basic_gas_viz Plots the concentration time curves of of CH4 and CO2 from
%enclosures
%   INPUT:
%       - ChamON_data:
%       - nchams: The number of chamber enclosures
%       - site_tag: Four letter code in all caps that designates where the
%       gases were collected
%       - ddmmmyyyy: Two-digit day, Three letter month, and four-digit year
%
%   OUTPUT:
%       - fig03: CH4 and CO2 concentration-time plots together in an array
%       - fig04: CH4 concentration-time plots together in an array
%       - fig05: CO2 concentration-time plots together in an array
%
%% Prelinary Visualizations

%-----Figure 03-----%
% Create a (2*dims_row) x (dims_col) array scatter plots of concentration
% vs time (rel sec) for both gases
fig03 = figure;
for idx = 1:nchams
    
%--CH4--%
    nexttile
    scatter(ChamON_data(:,2,idx),ChamON_data(:,3,idx));
    % This if-elseif-else block is used to appropriately designate the
    % numerical portion of each title in the array of plots. This block may
    % be removed if it is later determined to be unnecessary [24Feb2021]
            if idx     <= nchams
                trans   = 1;
                pnt     = idx;
            else
                trans   = 99;
                pnt     = 99;
            end
            % Title string (to be interated so the proper graph is
            % represented)
            title_str = sprintf(                                        ...
                       '%s %d.%d %s [CH_4] vs. Rel Time',               ...
                        site_tag, trans, pnt, ddmmmyyyy);
            title(title_str, 'FontSize', 8)
            % The maximum value on the x-axis will be 1200 (= 60 sec *
            % 20 min)
            xlim([0 1200])
            % Provide tick marks in 300 sec (5 min) intervals
            xticks(0:300:1200)
            % Add labels to axis
            xlabel('Relative Time (sec)', 'FontSize', 8)
            ylabel('[CH_4] (ppm)', 'FontSize', 8)
            % Apply gridding on both axis to see the data better
            grid on
            
%--CO2--%
    nexttile
    s_CO2 = scatter(ChamON_data(:,2,idx),ChamON_data(:,7,idx));
    s_CO2.MarkerEdgeColor = 'm';
    % This if-elseif-else block is used to appropriately designate the
    % numerical portion of each title in the array of plots
            if idx     <= nchams
                trans   = 1;
                pnt     = idx;
            else
                trans   = 99;
                pnt     = 99;
            end
            % Title string (to be interated so the proper graph is
            % represented)
            title_str = sprintf(                                        ...
                       '%s %d.%d %s [CO_2] vs. Rel Time',               ...
                        site_tag, trans, pnt, ddmmmyyyy);
            title(title_str, 'FontSize', 8)
            % The maximum value on the x-axis will be 1200 (= 60 sec *
            % 20 min)
            xlim([0 1200])
            % Provide tick marks in 300 sec (5 min) intervals
            xticks(0:300:1200)
            % Add labels to axis
            xlabel('Relative Time (sec)', 'FontSize', 8)
            ylabel('[CO_2] (ppm)', 'FontSize', 8)
            % Apply gridding on both axis to see the data better
            grid on        
                
end
% Save figure as .fig to working directory
fi       = sprintf("MATLAB_figs\\%s_%s_CH4_and_CO2_timeline_array.fig", ...
                    site_tag, ddmmmyyyy);
fig_file = working_dir+fi;
savefig(fig03, fig_file)

% Create a dims_row x dims_col array scatter plots of concentration vs time
% (rel sec)

% ---- Figure 04 ---- %
    % CH4
fig04 = figure;
for idx = 1:nchams

    nexttile
        scatter(ChamON_data(:,2,idx),ChamON_data(:,3,idx));
        % This if-elseif-else block is used to appropriately designate the
        % numerical portion of each title in the array of plots
                if idx     <= nchams
                    trans   = 1;
                    pnt     = idx;
                else
                    trans   = 99;
                    pnt     = 99;
                end
                % Title string (to be interated so the proper graph is
                % represented)
                title_str = sprintf(                                    ...
                       '%s %d.%d %s [CH_4] vs. Rel Time',               ...
                        site_tag, trans, pnt, ddmmmyyyy);
                title(title_str, 'FontSize', 8)
                % The maximum value on the x-axis will be 1200 (= 60 sec *
                % 20 min)
                xlim([0 1200])
                % Provide tick marks in 300 sec (5 min) intervals
                xticks(0:300:1200)
                % Add labels to axis
                xlabel('Relative Time (sec)', 'FontSize', 8)
                ylabel('[CH_4] (ppm)', 'FontSize', 8)
                % Apply gridding on both axis to see the data better
                grid on
end

% Save figure as .fig to working directory
fi       = sprintf("MATLAB_figs\\%s_%s_CH4_timeline_array.fig",         ...
                    site_tag, ddmmmyyyy);
fig_file = working_dir+fi;
savefig(fig04, fig_file)

% ---- Figure 05 ---- %
    % CO2
fig05 = figure;
for idx = 1:nchams
    
    nexttile
        s_CO2 = scatter(ChamON_data(:,2,idx),ChamON_data(:,7,idx));
        s_CO2.MarkerEdgeColor = 'm';
        % This if-elseif-else block is used to appropriately designate the
        % numerical portion of each title in the array of plots
                if idx     <= nchams
                    trans   = 1;
                    pnt     = idx;
                elseif idx  > dims_col
                    trans   = 2;
                    pnt     = idx-dims_col;
                else
                    trans   = 99;
                    pnt     = 99;
                end
                % Title string (to be interated so the proper graph is
                % represented)
                title_str = sprintf(                                    ...
                           '%s %d.%d %s [CO_2] vs. Rel Time',           ...
                            site_tag,trans,pnt, ddmmmyyyy); 
                title(title_str, 'FontSize', 8)
                % The maximum value on the x-axis will be 1200 (= 60 sec *
                % 20 min)
                xlim([0 1200])
                % Provide tick marks in 300 sec (5 min) intervals
                xticks(0:300:1200)
                % Add labels to axis
                xlabel('Relative Time (sec)', 'FontSize', 8)
                ylabel('[CO_2] (ppm)', 'FontSize', 8)
                % Apply gridding on both axis to see the data better
                grid on        
end

% Save figure as .fig to working directory
fi       = sprintf("MATLAB_figs\\%s_%s_CO2_timeline_array.fig",         ...
                    site_tag, ddmmmyyyy);
fig_file = working_dir+fi;
savefig(fig05, fig_file)

end


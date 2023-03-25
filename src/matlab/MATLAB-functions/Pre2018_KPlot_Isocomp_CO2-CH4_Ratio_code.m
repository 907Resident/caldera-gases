%% Data Visualization

%----Methane----%
% Plot regressions with 95% confidence intervals for the slope
fullscreen,
minimum = min(min(VCAC_11Jul2017_ISM_DATA(:,8,:)));
maximum = max(max(VCAC_11Jul2017_ISM_DATA(:,8,:)));
for i = 1:nchams
[y, yci]            = predict(CH4_lin_mdl{i}, VCAC_11Jul2017_ISM_DATA(:,1,i), 'Simultaneous', 'on');
coeff_slope         = table2array(CH4_lin_mdl{1,i}.Coefficients(2,1));
coeff_incp          = table2array(CH4_lin_mdl{1,i}.Coefficients(1,1));
r_squared           = CH4_lin_mdl{1,i}.Rsquared.Ordinary;
coeff_ci            = coefCI(CH4_lin_mdl{i});
regress_ttl_string  = sprintf('VCAC_11Jul2017 | Lin. Model | Trans 01 | Pt. %d', i);
fit_output          = sprintf(' [CH_4] = %2.3E * t + %2.3E \n Slope CI (%2.3f, %2.3f) \n Y-Intercept CI (%2.3f, %2.3f) \n R-squared = %1.4f', ...
                      coeff_slope, coeff_incp,                          ...
                      coeff_ci(1), coeff_ci(2),                         ...
                      coeff_ci(3), coeff_ci(4),                         ...
                      r_squared);
subplot(3,2,i)
    plot(VCAC_11Jul2017_ISM_DATA(:,1,i), y, 'r');
        hold on,
    plot(VCAC_11Jul2017_ISM_DATA(:,1,i), yci(:,1), '--r');
        hold on,
    plot(VCAC_11Jul2017_ISM_DATA(:,1,i), yci(:,2), '--r');
        hold on,
    plot(VCAC_11Jul2017_ISM_DATA(:,1,i), VCAC_11Jul2017_ISM_DATA(:,8,i), 'o');
        grid on
        title(regress_ttl_string, 'FontSize', 8)
        xlabel('Time (hr)', 'FontSize', 7)
        ylabel('[CH_4] (mg)', 'FontSize', 7)
        ylim([minimum maximum])
        xlim([0 0.25])
    if      i == 1
       annotation('textbox', [0.308 0.714 0.150 0.111], 'String', fit_output,  ...
                  'FitBoxToText','on', 'FontSize', 8, 'EdgeColor', 'None')
    elseif  i == 2
       annotation('textbox', [0.754 0.709 0.151 0.112], 'String', fit_output,  ...
                  'FitBoxToText','on', 'FontSize', 8, 'EdgeColor', 'None')
    elseif  i == 3
        annotation('textbox', [0.127 0.516 0.151 0.112], 'String', fit_output,  ...
                  'FitBoxToText','on', 'FontSize', 8, 'EdgeColor', 'None')
    end
end
clearvars coeff_slope coeff_incp r_squared coeff_ci regress_ttl_string fit_output maximum minimum

%----Carbon Dioxide----%
% Plot regressions with 95% confidence intervals for the slope
fullscreen,
minimum = min(min(VCAC_11Jul2017_ISM_DATA(:,9,:)));
maximum = max(max(VCAC_11Jul2017_ISM_DATA(:,9,:)));
for i = 1:nchams
[y, yci]            = predict(CO2_lin_mdl{i}, VCAC_11Jul2017_ISM_DATA(:,1,i), 'Simultaneous', 'on');
coeff_slope         = table2array(CO2_lin_mdl{1,i}.Coefficients(2,1));
coeff_incp          = table2array(CO2_lin_mdl{1,i}.Coefficients(1,1));
r_squared           = CO2_lin_mdl{1,i}.Rsquared.Ordinary;
coeff_ci            = coefCI(CO2_lin_mdl{i});
regress_ttl_string  = sprintf('VCAC_11Jul2017 | Lin. Model | Trans 01 | Pt. %d', i);
fit_output          = sprintf(' [CO_2] = %2.3E * t + %2.3E \n Slope CI (%2.3f, %2.3f) \n Y-Intercept CI (%2.3f, %2.3f) \n R-squared = %1.4f', ...
                      coeff_slope, coeff_incp,                          ...
                      coeff_ci(1), coeff_ci(2),                         ...
                      coeff_ci(3), coeff_ci(4),                         ...
                      r_squared);
subplot(3,2,i)
    plot(VCAC_11Jul2017_ISM_DATA(:,1,i), y, 'r');
        hold on,
    plot(VCAC_11Jul2017_ISM_DATA(:,1,i), yci(:,1), '--r');
        hold on,
    plot(VCAC_11Jul2017_ISM_DATA(:,1,i), yci(:,2), '--r');
        hold on,
    plot(VCAC_11Jul2017_ISM_DATA(:,1,i), VCAC_11Jul2017_ISM_DATA(:,9,i), 'mo');
        grid on
        title(regress_ttl_string, 'FontSize', 8)
        xlabel('Time (hr)', 'FontSize', 7)
        ylabel('[CH_4] (ppm)', 'FontSize', 7)
        ylim([minimum maximum])
        xlim([0 0.25])
    if      i == 1
       annotation('textbox', [0.308 0.714 0.150 0.111], 'String', fit_output,  ...
                  'FitBoxToText','on', 'FontSize', 8, 'EdgeColor', 'None')
    elseif  i == 2
       annotation('textbox', [0.754 0.709 0.151 0.112], 'String', fit_output,  ...
                  'FitBoxToText','on', 'FontSize', 8, 'EdgeColor', 'None')
    elseif  i == 3
        annotation('textbox', [0.127 0.516 0.151 0.112], 'String', fit_output,  ...
                  'FitBoxToText','on', 'FontSize', 8, 'EdgeColor', 'None')
    elseif i == 4
        annotation('textbox', [0.568 0.518 0.151 0.112], 'String', fit_output,  ...
                  'FitBoxToText','on', 'FontSize', 8, 'EdgeColor', 'None')
    elseif i == 5
        annotation('textbox', [0.324 0.117 0.151 0.112], 'String', fit_output,  ...
                  'FitBoxToText','on', 'FontSize', 8, 'EdgeColor', 'None')
    end
end
clearvars coeff_slope coeff_incp r_squared coeff_ci regress_ttl_string fit_output maximum minimum y yci

    %-----Exponential Model-----%
%     e_flux_time = cell([length(CH4_exp_flux), 1]);
%     T_exp_flux          = table;
%     for i = 1:nchams
%         flujo_tiempo    = datestr(VCAC_11Jul2017_ISM_DATA(490,1,i) + datenum(2017,03,06), 0);
%         e_flux_time{i}  = datetime(flujo_tiempo,'InputFormat','dd-MM-yyyy HH:mm:ss');
%     end
%
%     T_exp_flux.Time     = [e_flux_time{:,1}]';
%     T_exp_flux.Flux     = CH4_exp_flux(1,:)';
%
% fullscreen,
%     plot(T_exp_flux.Time, T_exp_flux.Flux, 'LineStyle', '--',           ...
%         'Marker', 'd');
%             title('CH_{4} Fluxes (Lin. Model) | BHCS')
%             ylabel('Flux (mg m^{-2} hr^{-1}')
%
% % Intercepts vs. Time
%
%     incp = zeros([nchams, 1]);
% for i = 1:nchams
%    incp(i) = CH4_KC_coeffs(1,2,i);
% end
%
% T_incp = table;
%     T_incp.Time = [e_flux_time{:,1}]';
%     T_incp.Incp = incp;
%         clearvars incp
%
% fullscreen,
%     plot(T_incp.Time, T_incp.Incp, 'LineStyle', '--', 'Marker', 'X', ...
%          'MarkerSize', 14);
%             title('Isotopic Composition of Methane | BHCS','Interpreter', 'tex')
%             ylabel('\delta^{13}C-CH_{4} (‰)','Interpreter', 'tex')
%             % Add the error at each point after the fact with textboxes

% CO2/CH4 Ratio per location
CO2_CH4_ratio = zeros([length(VCAC_11Jul2017_ISM_DATA), 1]);
for i = 1:nchams
    CO2_CH4_ratio(:,i) =  VCAC_11Jul2017_ISM_DATA(:,6,i)./VCAC_11Jul2017_ISM_DATA(:,2,i);
end

%% CO2/CH4 Ratio v. Time
fullscreen,
minimum = min(min(CO2_CH4_ratio));
maximum = max(max(CO2_CH4_ratio));
for i = 1:nchams
    title_string = sprintf('VCAC_11Jul2017 | CO_2:CH_4 Ratio vs. Time | Trans 01 | Pt %d', i);
    subplot(3,2,i)
        scatter(VCAC_11Jul2017_ISM_DATA(:,1,i), CO2_CH4_ratio(:,i));
        title(title_string)
        ylabel('[CO_2]:[CH_4] Ratio')
        xlabel('Time (hr)')
        ylim([minimum maximum]);
        xlim([0 0.25])
end
clearvars title_string maximum minimum

% CO2:CH4 Flux Ratio (mg m-2 hr-1)
    % Linear Model
Flx_Ratio = CO2_lin_flux ./ CH4_lin_flux;
    % Visualize
    fullscreen,
    subplot(1,2,1)
    b1 = bar(Flx_Ratio, 'FaceColor', [0.3010 0.7450 0.9330],            ...
            'EdgeColor', 'k');
            title('CO_2-CH_4 Flux Ratio | VCAC_11Jul2017 | Trans 01')
            ylabel('F_{CO2_2}:F_{CH_4} (mg m^{-2} hr^{-1})')
            xticklabels({'Point 01', 'Point 02', 'Point 03', 'Point 04',...
                         'Point 05'});
    subplot(1,2,2)
    b2 = bar(log10(Flx_Ratio), 'FaceColor', [0.3010 0.7450 0.9330],            ...
            'EdgeColor', 'k');
            ylabel('log_{10}(F_{CO2_2}:F_{CH_4}) (mg m^{-2} hr^{-1})')
            xticklabels({'Point 01', 'Point 02', 'Point 03', 'Point 04',...
                         'Point 05'});


%% Compare carbon isotopes for both gases
VCAC_11Jul2017_isocomp    = zeros([1 nchams]);
for i = 1:nchams
    % Set up the vectors for that exclude the NaNs
    x = ~isnan(VCAC_11Jul2017_ISM_DATA(:,4,i));
        x = find(x);
    xData = VCAC_11Jul2017_ISM_DATA(1:x(end),4,i);
    y = ~isnan(VCAC_11Jul2017_ISM_DATA(:,7,i));
        y = find(y);
    yData = VCAC_11Jul2017_ISM_DATA(1:y(end),7,i);

    VCAC_11Jul2017_isocomp(i) = isocomp(xData,yData);
    % Create a vector of RGB values for time gradient inset
    ce = length(xData);
    cf = jet(ce);
    fullscreen,
    h  = scatter(xData,yData, 14, cf, 'filled');
         if i <= 4
            scatt_title_string =                                        ...
                sprintf('^{13}C-CO_2 vs. ^{13}C-CH_4 | Trans 01 | Pt. %d', i);
         else
             scatt_title_string =                                        ...
                 sprintf('^{13}C-CO_2 vs. ^{13}C-CH_4 | Trans 02 | Pt. %d', i-4);
         end
         title(scatt_title_string)
         colormap jet
         axis([min(VCAC_11Jul2017_ISM_DATA(:,4,i))-2 max(VCAC_11Jul2017_ISM_DATA(:,4,i))+2,      ...
               min(VCAC_11Jul2017_ISM_DATA(:,7,i))-2 max(VCAC_11Jul2017_ISM_DATA(:,7,i))+2]);
         set(gca, 'Color', [0.05 0.15 0.15])
         grid on
end

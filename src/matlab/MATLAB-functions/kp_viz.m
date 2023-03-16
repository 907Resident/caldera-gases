function [fig09,fig10] = kp_viz(ChamON_data, nchams,                    ...
                                KP_Fits,                                ...
                                KP_coeffs, KP_coeffs_ci,                ...
                                site_tag, ddmmmyyyy,                    ...
                                matlab_fig_dir)
%KP_VIZ Plots Keeling Plot intercepts and an array of full Keeling Plots as
%regressions
%   Detailed explanation goes here
%% Visualize the Keeling Plot Intercepts
    % Add the error bars for appropriate sites
fig09 = figure;
    %  CH4
        nexttile
            dmy_Y_CH4     = KP_coeffs(:,2,:,1);
            Y             = reshape(dmy_Y_CH4,[1 nchams]);
            neg_err_CH4   = reshape(KP_coeffs_ci(1,2,:,1),[1 nchams]);
            pos_err_CH4   = reshape(KP_coeffs_ci(2,2,:,1),[1 nchams]);
            bar(Y)
            hold on
            errorbar(1:nchams,Y,abs(Y - neg_err_CH4),                   ...
                                abs(Y - pos_err_CH4),                   ...
                    'LineStyle', 'none', 'Color', 'k')
            grid on
            title_str = sprintf('%s %s CH_4 Keeling Plot | \\delta^{13}C-CH_4 Signatures',...
                                 site_tag, ddmmmyyyy);
            title(title_str, 'FontSize', 12,'Interpreter', 'tex') 
            xlabel('Location')
            ylabel('\delta^{13}C-CH_4 (‰)')
                     
    % CO2
        nexttile
            dmy_Y_CO2     = KP_coeffs(:,2,:,2);
            Y             = reshape(dmy_Y_CO2,[1 nchams]);
            neg_err_CO2   = reshape(KP_coeffs_ci(1,2,:,2), [1 nchams]);
            pos_err_CO2   = reshape(KP_coeffs_ci(2,2,:,2), [1 nchams]);
            bar(Y, 'FaceColor', 'm')
            hold on
            errorbar(1:nchams,Y,abs(Y - neg_err_CO2),                   ...
                                abs(Y - pos_err_CO2),                   ...
                    'LineStyle', 'none', 'Color', 'k')
            grid on
            title_str = sprintf('%s %s CO_2 Keeling Plot | \\delta^{13}C-CO_2 Signatures',...
                                 site_tag, ddmmmyyyy);
            title(title_str, 'FontSize', 12,'Interpreter', 'tex')            
            xlabel('Location')
            ylabel('\delta^{13}C-CO_2 (‰)','Interpreter', 'tex')
 
% Save figure as .fig to working directory
fi       = sprintf("%s_%s_CH4_and_CO2_d13C_KP_sig_bar.fig",...
                    site_tag, ddmmmyyyy);
fig_file = fullfile(matlab_fig_dir, fi);
savefig(fig09, fig_file)

% Visualize Keeling Plots in an array of regressions
fig10 = figure;
for idx = 1:nchams
%--CH4--%
    nexttile
        plot(KP_Fits{idx,1},                                            ...
                  1./ChamON_data(:,3,idx),ChamON_data(:,5,idx),         ...
                  'predfunc')
        fxn        = sprintf('δ^{13}C-CH_4 = %2.2f*[CH_4]^{-1} + %2.2f',...
                              KP_coeffs(1,1,idx,1),KP_coeffs(1,2,idx,1));
        str        = sprintf('R^{2} = %1.3f', KP_Fits{idx,2});
        annotation('textbox',[0 0.8 0.5 0.2],'String',fxn,              ...
                   'FitBoxtoText','on', 'FontSize', 7,                  ...
                   'FontWeight', 'bold',                                ...
                   'BackgroundColor', [0.97 0.97 0.97], 'FaceAlpha',0.5,...
                   'EdgeColor', 'none')
        annotation('textbox',[0 0.5 0.3 0.2],'String',str,              ...
                   'FitBoxtoText','on' , 'FontSize', 7)
        xlabel('[CH_4]^{-1} (ppm^{-1})','FontSize', 7)
        ylabel('\delta^{13}C-CH_4 (‰)','Interpreter', 'tex',           ...
               'FontSize', 7)
        grid on
    % This if-elseif-else block is used to appropriately designate the
    % numerical portion of each title in the array of plots
                if  idx    <= nchams
                    trans   = 1;
                    pnt     = idx;
                else
                    trans   = 99;
                    pnt     = 99;
                end
            % Title string (to be interated so the proper graph is
            % represented)
            title_str = sprintf(                                        ...
                       '%d.%d CH_4 Keeling Plot',                       ...
                        trans,pnt);
            title(title_str, 'FontSize', 8)
%--CO2--%
    nexttile
        p_CO2 = plot(KP_Fits{idx,3},                                    ...
                1./ChamON_data(:,7,idx),ChamON_data(:,8,idx),'predfunc');
        p_CO2(1).Color = 'm';
        p_CO2(2).Color = 'k';
        p_CO2(3).Color = 'k';
        fxn            = sprintf('δ^{13}C-CO_2 = %2.2f*[CO_2]^{-1} + %2.2f',...
                                 KP_coeffs(1,1,idx,2),KP_coeffs(1,2,idx,2));                    
        str            = sprintf('R^{2} = %1.3f', KP_Fits{idx,4});
        annotation('textbox',[0 0.2 0.5 0.2],'String',fxn,              ...
                   'FitBoxtoText','on', 'FontSize', 7,                  ...
                   'FontWeight', 'bold',                                ...
                   'BackgroundColor', [0.97 0.97 0.97], 'FaceAlpha',0.5,...
                   'EdgeColor', 'none')
        annotation('textbox',[0 0 0.3 0.2],'String',str,                ...
                   'FitBoxtoText', 'on', 'FontSize', 7)
        xlabel('[CO_2]^{-1} (ppm^{-1})','FontSize',7)
        ylabel('\delta^{13}C-CO_2 (‰)','FontSize',7)
        grid on
    % This if-elseif-else block is used to appropriately designate the
    % numerical portion of each title in the array of plots
                if  idx    <= nchams
                    trans   = 1;
                    pnt     = idx;
                else
                    trans   = 99;
                    pnt     = 99;
                end
            % Title string (to be interated so the proper graph is
            % represented)
            title_str = sprintf('%s %d.%d %s CO_2 Keeling Plot',        ...
                                 site_tag, trans, pnt, ddmmmyyyy);
            title(title_str, 'FontSize', 8)
end

% Save figure as .fig to working directory
fi       = sprintf("%s_%s_CH4_and_CO2_d13C_KP_sig_bar.fig",...
                    site_tag, ddmmmyyyy);
fig_file = fullfile(matlab_fig_dir, fi);
savefig(fig10, fig_file)
end


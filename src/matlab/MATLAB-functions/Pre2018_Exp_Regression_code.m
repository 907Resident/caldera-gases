%% Old Exponential Regression Code

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
%    incp(i) = KC_coeffs(1,2,i);
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

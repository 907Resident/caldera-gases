function [scat_plot, sp, iso_array] = isocomp_array( ...
    iCH4,iCO2,nchams, varargin)
%ISOCOMP_ARRAY(iCH4,iCO2,dims)
%  iCH4:        isotopic composition of methane
%  iCO2:        isotopic composition of carbon dioxide
%  nchams:      number of chamber measurements in the array (cannot exceed
%               the product of dims_rows and dims_cols
%  varargin:    Enter if you want the plot zoomed in - "Zoom", "yes"
%% Plot Isocomps e.g., Whiticar (1999)

% Conditional block for reshaping a > 2D array of isotopic data
iso_dims_CH4    = size(iCH4);
iso_dims_CO2    = size(iCO2);

if length(iso_dims_CH4) > 2
   iCH4         = reshape(iCH4,[length(iCH4),nchams]);
elseif length(iso_dims_CO2) > 2
   iCO2         = reshape(iCO2,[length(iCO2),nchams]);
end

% Preallocate scatter plots (scat_plot) and subplots (sp)
scat_plot       = cell([nchams 1]);
sp              = cell([nchams 1]);

% Create figure
iso_array       = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
tiles           = tiledlayout('flow');
for plc         = 1:nchams

% Create scatter
    % Create vector for color gradient to represent time
    ce          = find(isnan(iCH4(:,plc)),1,'first');
    sz          = size(ce);
        if sz(1)== 0
           ce = find(~isnan(iCO2(:,1,plc)), 1, 'last');
        end
    cf          = jet(ce);
    nexttile(tiles)
    scat_plot   = scatter(iCH4(1:ce,plc),iCO2(1:ce,plc), 14,            ...
                          cf, 'filled');
        colormap jet
        % Create colorbar for the last plot in the array
        if plc == nchams
            % Create the colorbar which will enable the interpretation of
            % the colored markers on the plot
            cbar    = colorbar('eastoutside',                           ...
                               'Ticks', [0 0.5 1.0],                    ...
                               'TickLabels',{'Beginning','Middle','End'});
            cbar.Label.String ='Relative Time of Chamber Measurement';
            % Rotatate the label on the colorbar so that it is legigble
            set(get(cbar,'Label'),'Rotation', 270.0)
        end
        % Set the limits of the x and y axis
        if strcmp(varargin{2}, "yes") == 1
            xlim_lwr = 5 .* floor(nanmin(iCH4(:,plc)) ./ 5);
            xlim_upr = 5 .* ceil(nanmax(iCH4(:,plc))  ./ 5);
            ylim_lwr = 5 .* floor(nanmin(iCO2(:,plc)) ./ 5);
            ylim_upr = 5 .* ceil(nanmax(iCO2(:,plc))  ./ 5);        
            xlim([xlim_lwr xlim_upr])
            ylim([ylim_lwr ylim_upr])
        else 
            xlim([-100 -20])
            ylim([-40   20])
        end
        % Label the x and y axis
        xlabel('\delta^{13}C-CH_{4} (‰)','Interpreter','tex',          ...
               'FontSize', 11)
        ylabel('\delta^{13}C-CO_{2} (‰)','Interpreter','tex',          ...
               'FontSize', 11)
        % Create iterative title string and title the graphs
        ttl_str = sprintf('\x03B4^{13}C-CO_2 vs \x03B4^{13}C-CH_4 Pt# %d', ...
                           plc);
        title(ttl_str,'FontSize',10, 'Interpreter', 'tex')
        
        hold on
        
% Create and plot fractionation lines
    % These are based off Eq. (7) in Whiticar (1999)
        % Fractionation Line #1 episolon = 90
        ex = -100:10:-70;
        ey = -010:10:020;
        plot(ex,ey,'-k','LineWidth',0.8)
        
        % Fractionation Line #2 epislon = 80
        ex = -100:10:-60;
        ey = -020:10:020;
        plot(ex,ey,'-k','LineWidth',0.8)
        
        % Fractionation Line #2 epislon = 70
        ex = -100:10:-50;
        ey = -030:10:020;
        plot(ex,ey,'-k','LineWidth',0.8)
        
        % Fractionation Line #4 epislon = 60
        ex = -100:10:-40;
        ey = -040:10:020;
        plot(ex,ey,'-k','LineWidth',0.8)
        
        % Fractionation Line #5 epislon = 50
        ex = -100:10:-30;
        ey = -050:10:020;
        plot(ex,ey,'-k','LineWidth',0.8)
        
        % Fractionation Line #6 epislon = 40
        ex = -090:10:-20;
        ey = -050:10:020;
        plot(ex,ey,'-k','LineWidth',0.8)
        
        % Fractionation Line #7 epislon = 30
        ex = -080:10:-20;
        ey = -050:10:010;
        plot(ex,ey,'-k','LineWidth',0.8)
        
        % Fractionation Line #8 epislon = 20
        ex = -070:10:-20;
        ey = -050:10:000;
        plot(ex,ey,'-k','LineWidth',0.8)
        
        % Fractionation Line #9 epislon = 5
        ex = -055:10:-15;
        ey = -050:10:-10;
        plot(ex,ey,'-k','LineWidth',0.8)

        % Set the remaining axes properties
        set(gca,'XGrid','on','XMinorTick','on','YGrid','on',        ...
                    'YMinorTick','on');

        hold on
        
% Create a diamond to indicate where the atmosphere is located for both
% isotopes
        scatter(-47,-9,80,'diamond', 'MarkerFaceColor', 'k')        
        
end
        

end
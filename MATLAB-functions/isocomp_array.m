function [scat_plot, sp] = isocomp_array(iCH4,iCO2,nchams)
%ISOCOMP_ARRAY(iCH4,iCO2,dims)
%  iCH4:        isotopic composition of methane
%  iCO2:        isotopic composition of carbon dioxide
%  S1:          scatter s
%  C1:          scatter c
%  dims_row:    number of rows in the array (must be greater than zero)
%  dims_col:    number of columns in the array (must be greater than zero)
%  nchams:      number of chamber measurements in the array (cannot exceed
%               the product of dims_rows and dims_cols
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
%iso_array       = figure;
fullscreen
tiles           = tiledlayout('flow');
for plc         = 1:nchams
% Create subplot
%sp{plc}         = subplot(dims_row,dims_col,plc,'Parent',iso_array);

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
        % Create the colorbar which will enable the interpretation of the
        % colored markers on the plot
        cbar    = colorbar('eastoutside',                               ...
                           'Ticks', [0 0.5 1.0],                        ...
                           'TickLabels',{'Beginning','Middle','End'});
        cbar.Label.String ='Relative Time of Chamber Measurement';
        % Rotatate the label on the colorbar so that it is legigble
        set(get(cbar,'Label'),'Rotation', 270.0)
        % Set the limits of the x and y axis
        xlim([-100 -20]), ylim([-50 20])
        % Label the x and y axis
        xlabel('\delta^{13}C-CH_{4} (‰)','Interpreter','tex', ...
               'FontSize', 11)
        ylabel('\delta^{13}C-CO_{2} (‰)','Interpreter','tex', ...
               'FontSize', 11)
        % Create iterative title string and title the graphs
        ttl_str = sprintf('\x03B4^{13}C-CO_2 vs \x03B4^{13}C-CH_4 Pt# %d',...
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


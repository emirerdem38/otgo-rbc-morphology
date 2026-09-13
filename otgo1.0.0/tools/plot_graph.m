function [] = plot_graph(image_title,div,data,analysis_type,saveto,varargin)

% PLOTGRAPH plots the forces and torques for an analysis, if an interval is
% provided, it calculates the slope for the specified data

if ~isempty(varargin)
    interval = varargin{1}; 
    order = varargin{2};
    try
        slope_order = varargin{3};
    catch
    end

    data_points_x = div((div >= interval(1) & div <= interval(2)));
    data_points_y = data(order, (div >= interval(1) & div <= interval(2)));
    coefficients = polyfit(data_points_x, data_points_y, 1); % Perform linear fit
    slope = coefficients(1);
end

switch analysis_type
    case 'X_disp'
        label = 'x';
    case 'Y_disp'
        label = 'y';
    case 'Z_disp'
        label = 'z';
    case 'X_rot'
        label = 'x';
    case 'Y_rot'
        label = 'y';
    case 'Z_rot'
        label = 'z';
end

if isempty(varargin)

    % Plot Forces
    figure; 
    hold on 
    box on;
    %title(sprintf('%s %s', '(Forces)', image_title))
    plot(div,data(1,:),'LineWidth',2)
    plot(div,data(2,:),'LineWidth',2)
    plot(div,data(3,:),'LineWidth',2)
    if ismember(analysis_type, {'X_disp','Y_disp','Z_disp'})
        %xlabel(sprintf('%s coordinates (m)', upper(label)),'FontSize',18)
        xlabel('\textbf{$\mu m$}','Interpreter','latex','FontSize',24);
    else
        xlabel(sprintf('Degrees of rotation along %s-axis (rad)',label),'FontSize',24)
    end
    %ylabel('Forces (pN)','FontSize',18)
    ylabel('\textbf{F ($pN$)}','Interpreter','latex','FontSize',24);
    line([min(div) max(div)], [0 0], 'Color', 'k', 'LineStyle', '--','HandleVisibility', 'off')
    legend('Fx','Fy','Fz','Location','best', 'Box', 'off')
    % set(gca,'linewidth',2)
    set(gca,'linewidth',1, 'fontsize' , 22)
    try
        saveas(gcf, sprintf('%s%s%s', saveto, 'Forces', '.jpg'))
    catch
    end
    
    % Plot Torques 
    figure;
    hold on
    box on;
    %title(sprintf('%s %s', '(Torques)', image_title))
    plot(div,data(4,:),'LineWidth',2)
    plot(div,data(5,:),'LineWidth',2)
    plot(div,data(6,:),'LineWidth',2)
    if ismember(analysis_type, {'X_disp','Y_disp','Z_disp'})
        xlabel(sprintf('%s coordinates (m)', upper(label)),'FontSize',24)
    else
        xlabel(sprintf('Degrees of rotation along %s-axis (rad)',label),'FontSize',24)
    end
    ylabel('Torques (Nm)','FontSize',24)
    line([min(div) max(div)], [0 0], 'Color', 'k', 'LineStyle', '--','HandleVisibility', 'off')
    legend('Tx','Ty','Tz','Location','best', 'Box', 'off')
    set(gca,'linewidth',2)
    try
        saveas(gcf, sprintf('%s%s%s', saveto, 'Torques', '.jpg'))
    catch
    end

else

    % Plot single entity with slope
    figure;
    hold on
    box on;
    outside_x = div((div < interval(1) | div > interval(2)));
    outside_y = data(order, (div < interval(1) | div > interval(2)));

    plot(outside_x,outside_y,'LineWidth',2,'LineStyle','--','Color','b')
    
    plot(data_points_x, polyval(coefficients, data_points_x), 'r', 'LineWidth', 3) % Add linear fit line
    if order <= 3
        %title(sprintf('(%s Force) %s (Linear Fit)', upper(label), image_title))
        xlabel(sprintf('%s coordinates (m)', upper(label)),'FontSize',18)
        ylabel(sprintf('%s Force (N)',upper(label)),'FontSize',18)
        legend(sprintf('F%s',label),'Linear Fit','Location','best', 'Box', 'off')
    else
        %title(sprintf('(%s Torque) %s (Linear Fit)', upper(label), image_title))
        xlabel(sprintf('Degrees of rotation along %s-axis (rad)',label),'FontSize',18)
        ylabel(sprintf('%s Torque (Nm)',upper(label)),'FontSize',18)
        legend(sprintf('T%s',label),'Linear Fit','Location','best', 'Box', 'off')  
    end
    
    set(gca,'linewidth',2)
    text_x = (interval(1) + interval(2)) / 2;
    text_y = polyval(coefficients, text_x);
    disp(slope)
    text(text_x * (interval(2) - interval(1)), text_y + 0.2 * (max(data_points_y) - min(data_points_y)), sprintf('Slope:\n %.2fx10^{-18} Nm/rad', slope*1e+18), 'FontSize', 15, 'Color', 'r', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
    
    try
        if order <= 3
            try
                saveas(gcf, sprintf('%s(%s Force)(Linear Fit)(%d).jpg', saveto, upper(label),slope_order))
            catch
                saveas(gcf, sprintf('%s(%s Force)(Linear Fit).jpg', saveto, upper(label)))
            end
        else
            try
                saveas(gcf, sprintf('%s(%s Torque)(Linear Fit)(%d).jpg', saveto, upper(label),slope_order))
            catch
                saveas(gcf, sprintf('%s(%s Torque)(Linear Fit).jpg', saveto, upper(label)))
            end
        end
    catch
    end
end



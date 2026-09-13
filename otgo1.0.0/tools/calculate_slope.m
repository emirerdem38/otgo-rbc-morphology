function slope = calculate_slope(forces_torques,div,order,int_start,int_end)

% CALCULATE_SLOPE calculates the slope when the data of forces and torques
% are provided for a certain analysis.
%
% DIV is the total interval for which the analysis is made.
% forces_torques is the 6 by N matrix denoting all the forces and torques
% along 3 directions of N instances. START and END  denotes the interval
% for which the fit will be made. 

% switch data
%     case 'FX'
%         order = 1;
%     case 'FY'
%         order = 2;
%     case 'FZ'
%         order = 3;
%     case 'TX'
%         order = 4;
%     case 'TY'
%         order = 5;
%     case 'TZ'
%         order = 6;
% end

data_points_x = div((div >= int_start & div <= int_end));
data_points_y = forces_torques(order, (div >= int_start & div <= int_end));

coefficients = polyfit(data_points_x, data_points_y, 1); % Perform linear fit
slope = coefficients(1);

end


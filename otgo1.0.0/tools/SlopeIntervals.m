function intervals = SlopeIntervals(X,Y)

% Will return the possible intervals that a line can be fitted

enough_for_fit = numel(X)/10; % How many indexes is seen enough for a fit
filter_crossing_points = "yes";
crossing_points = ZeroCross(X,Y,filter_crossing_points);
slope_indexes = [];
exists = 0;
intervals = [];

if numel(crossing_points) ~= 0
   
    for i = crossing_points

        until_neg = i - StraightUntil(X,Y,i,-1);
        until_pos = StraightUntil(X,Y,i,1) - i;
    
        for j = slope_indexes

            if abs(i-j) < (enough_for_fit-1)
                exists = 1;
            end
        end

        if until_neg + until_pos + 1 >= enough_for_fit && exists == 0
            slope_indexes = [slope_indexes i];
            intervals = [intervals; X(i-until_neg) X(i + until_pos) i];
        end
        exists = 0;

    end
    
end

end


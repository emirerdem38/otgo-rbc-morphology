function crossing_points = ZeroCross(X,Y,filter,varargin)

if ~isempty(varargin)
    cross = NaN;
    interval = varargin{1}; 
end

crossing_points = [];
for i = 2:numel(Y)
    % Check if the current and previous points have opposite signs
    if sign(Y(i)) ~= sign(Y(i-1))
        % Use linear interpolation to approximate the crossing point
        crossing_point = interp1([Y(i-1), Y(i)], [X(i-1), X(i)], 0);
        crossing_points = [crossing_points, crossing_point];
    end
end

if filter == "yes" && numel(crossing_points) ~= 0

    idx = zeros(1,numel(crossing_points));

    for i = 1:numel(crossing_points)
        [~, idx(i)] = min(abs(X - crossing_points(i)));
    end

    thresholdDistance = numel(X)/20;
 
    selectedIndexes = idx(1); 
    lastSelectedIndex = idx(1);
    for i = 2:length(idx)
        % Check if the current index is outside the threshold distance from the last selected index
        if abs(idx(i) - lastSelectedIndex) >= thresholdDistance
            % If it is, add it to the selected indexes
            selectedIndexes = [selectedIndexes, idx(i)];
            lastSelectedIndex = idx(i);  % Update the last selected index
        end
    end

    crossing_points = selectedIndexes;
end

if ~isempty(varargin) && numel(crossing_points) ~= 0
    for i = 1:numel(crossing_points)
        if X(crossing_points(i)) > interval(1) && X(crossing_points(i)) < interval(2)
            cross = crossing_points(i);
        end
    end
    crossing_points = cross;
end

end


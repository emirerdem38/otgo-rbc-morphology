function diff = compare(cell,p_set)

mag = cell.mag;

differences = zeros(numel(p_set),1);
for i = 1:numel(p_set)
    dists = zeros(size(cell.cart.x,1),1);
    for j = 1:size(cell.cart.x,1)
        dists(j) = ( (cell.cart.x(j) - p_set.X(i)*10^(-mag))^2 + (cell.cart.y(j) - p_set.Y(i)*10^(-mag))^2 + (cell.cart.z(j) - p_set.Z(i)*10^(-mag))^2 )^(1/2);
    end
    differences(i) = min(dists);
end

diff = mean(differences);

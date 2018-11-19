function distance = get_distance_between_points(i, j, grid)

Infty = 1e+10;
k = get_point_index(i);
l = get_point_index(j);

distance = (norm(grid(k, :) - grid(l, :)))^2 * (rem(i, 3) == rem(j, 3)) + Infty * (rem(i, 3) ~= rem(j, 3));

end


function k = get_point_index(i)

k = ceil(i / 3);

end
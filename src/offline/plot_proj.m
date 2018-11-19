function plot_proj( restored_beta, grid, time, proj, moment )
% proj = 'xy', 'xz' or 'yz'
% moment = 1 or 2 or 3

q = restored_beta(:, time);
w = q(moment:3:end);

un_x = unique(grid(:, 1));
un_y = unique(grid(:, 2));
un_z = unique(grid(:, 3));

clims = [-600, 600];

if (strcmp(proj, 'xy'))
    space_beta = zeros(numel(un_x), numel(un_y));
    for i = 1 : numel(un_x)
        for j = 1 : numel(un_y)
            space_beta(i, j) = sum(w(grid(:, 1) == un_x(i) & grid(:, 2) == un_y(j)));
        end
    end
    imagesc(space_beta, clims)
    colorbar
    
elseif (strcmp(proj, 'xz'))
    space_beta = zeros(numel(un_x), numel(un_z));
    for i = 1 : numel(un_x)
        for j = 1 : numel(un_z)
            space_beta(i, j) = sum(w(grid(:, 1) == un_x(i) & grid(:, 3) == un_z(j)));
        end
    end
    imagesc(space_beta, clims)
    colorbar
    
elseif (strcmp(proj, 'yz'))
    space_beta = zeros(numel(un_y), numel(un_z));
    for i = 1 : numel(un_y)
        for j = 1 : numel(un_z)
            space_beta(i, j) = sum(w(grid(:, 2) == un_y(i) & grid(:, 3) == un_z(j)));
        end
    end
    imagesc(space_beta, clims)
    colorbar
    
end


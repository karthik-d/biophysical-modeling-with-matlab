%% (A) Uniform random points.
N = 50;
randomPosns = rand(2, N);

[V_random, C_random, vAll_random, cAll_random, xbox_random, ybox_random] = getVoronoiDiagram( ... 
    randomPosns(1, :), randomPosns(2, :), 1);

% save generated data and plot.
saveas(gcf, 'outputs/voronoi_random.fig'); clf;
writematrix(V_random, 'outputs/vertexPosns_random.dat');
writecell(C_random, 'outputs/cellVertices_random.dat');


%% (B-1): Disordered, jammed disk packing.
diskwallPosns = readmatrix('data/disks_walls.dat');

[V_disk, C_disk, vAll_disk, cAll_disk, xbox_disk, ybox_disk] = getVoronoiDiagram( ...
    diskwallPosns(:, 2), diskwallPosns(:, 3), 1);
% save plot.
saveas(gcf, 'outputs/voronoi_diskwall.fig'); clf;


%% (B-2): Lloyd's algorithm.
lloydPosns = readmatrix('data/lloyds.dat');

[V_lloyd, C_lloyd, vAll_lloyd, cAll_lloyd, xbox_lloyd, ybox_lloyd] = getVoronoiDiagram( ...
    lloydPosns(:, 2), lloydPosns(:, 3), 1);
% save plot.
saveas(gcf, 'outputs/voronoi_lloyd.fig'); clf;


%% (B-3): Stomata positions.
stomataPosns = readmatrix('data/stomata.dat');

[V_stomata, C_stomata, vAll_stomata, cAll_stomata, xbox_stomata, ybox_stomata] = getVoronoiDiagram( ...
    stomataPosns(:, 2), stomataPosns(:, 3), 1);
% save plot.
saveas(gcf, 'outputs/voronoi_stomata.fig'); clf;


% (C) Lewis law.
areas = plotLewisLaw(V_stomata, C_stomata, xbox_stomata, ybox_stomata);
print(areas);


function [cellAreas] = plotLewisLaw(vertexPosns, cellNeighbors, cellsX, cellsY)

    num_cells = numel(cellNeighbors);
    n_distribution = cellfun(@numel, cellNeighbors);
    n_bar = sum(n_distribution)/num_cells;
    
    % Accumulate cell areas.
    cellAreas = zeros(num_cells, 1);
    for i=1:num_cells
        vertices = vertexPosns(cellNeighbors{i}, :);
        
        totalArea = 0;
        for j=1:(size(vertices, 1)-1)
            totalArea = totalArea + ...
                triangularArea([cellsX(j); cellsY(j)], vertices(j, :), vertices(j+1, :));
        end
        cellAreas(i) = totalArea;

    end
end
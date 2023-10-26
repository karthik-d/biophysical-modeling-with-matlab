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


% Record cell perimeters.
nNaught = 2;

[nDistribution_random, nBar_random, cellAreas_random, cellPeris_random, cellsPerN_random, cellAreasPerN_random] = ...
	computeCellParams(V_random, C_random, xbox_random, ybox_random);

[nDistribution_disk, nBar_disk, cellAreas_disk, cellPeris_disk, cellsPerN_disk, cellAreasPerN_disk] = ...
	computeCellParams(V_disk, C_disk, xbox_disk, ybox_disk);

[nDistribution_lloyd, nBar_lloyd, cellAreas_lloyd, cellPeris_lloyd, cellsPerN_lloyd, cellAreasPerN_lloyd] = ...
	computeCellParams(V_lloyd, C_lloyd, xbox_lloyd, ybox_lloyd);

[nDistribution_stomata, nBar_stomata, cellAreas_stomata, cellPeris_stomata, cellsPerN_stomata, cellAreasPerN_stomata] = ...
	computeCellParams(V_stomata, C_stomata, xbox_stomata, ybox_stomata);


% (C) Lewis law.
cellAreasAvg_random = cellAreasPerN_random ./ cellsPerN_random;
cellAreasAvg_disk = cellAreasPerN_disk ./ cellsPerN_disk;
cellAreasAvg_lloyd = cellAreasPerN_lloyd ./ cellsPerN_lloyd;
cellAreasAvg_stomata = cellAreasPerN_stomata ./ cellsPerN_stomata;

figure(2); hold on;
plot(cellAreasAvg_random ./ cellAreasAvg_random(nBar_random), ...
	((1:length(cellsPerN_random))-nNaught) ./ (nBar_random-nNaught));
plot(cellAreasAvg_disk ./ cellAreasAvg_disk(nBar_disk), ...
	((1:length(cellsPerN_disk))-nNaught) ./ (nBar_disk-nNaught));
plot(cellAreasAvg_lloyd ./ cellAreasAvg_lloyd(nBar_lloyd), ...
	((1:length(cellsPerN_lloyd))-nNaught) ./ (nBar_lloyd-nNaught));
plot(cellAreasAvg_stomata ./ cellAreasAvg_stomata(nBar_stomata), ...
	((1:length(cellsPerN_stomata))-nNaught) ./ (nBar_stomata-nNaught));

disp(((1:length(cellsPerN_random))-nNaught) ./ (nBar_random-nNaught));
disp(cellAreasAvg_random ./ cellAreasAvg_random(nBar_random));


% Function to compute cell parameters.
function [nDistribution, nBar, cellAreas, cellPerimeters, cellsPerN, cellAreasPerN] = ...
    computeCellParams(vertexPosns, cellNeighbors, cellsX, cellsY)

    num_cells = numel(cellNeighbors);
    nDistribution = cellfun(@numel, cellNeighbors);
    nBar = round(sum(nDistribution)/num_cells);
    
    % Accumulate cell areas.
    cellAreasPerN = zeros(max(nDistribution), 1);
	cellsPerN = zeros(max(nDistribution), 1);
    cellAreas = zeros(num_cells, 1);
    cellPerimeters = zeros(num_cells, 1);
    for i=1:num_cells
        vertices = vertexPosns(cellNeighbors{i}, :);
        
        totalArea = 0;
        totalPerimeter = 0;
        for j=1:(size(vertices, 1)-1)
            totalPerimeter = totalPerimeter + ...
                pdist2(vertices(j, :), vertices(j+1, :));
            totalArea = totalArea + ...
                triangularArea([cellsX(i) cellsY(i)], vertices(j, :), vertices(j+1, :));
        end
        
		cellsPerN(nDistribution(i)) = cellsPerN(nDistribution(i)) + 1;
		cellAreasPerN(nDistribution(i)) = cellAreasPerN(nDistribution(i)) + totalArea;
        cellAreas(i) = totalArea;
        cellPerimeters(i) = totalPerimeter;
	end
end
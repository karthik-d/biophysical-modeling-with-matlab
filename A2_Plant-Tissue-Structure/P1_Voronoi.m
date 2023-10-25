N = 50;

randomPosns = rand(2, N);

figure(1);
[V, C, vAll, cAll, xbox, ybox] = getVoronoiDiagram(randomPosns(1, :), randomPosns(2, :), 1);

% save generated data and plot.
saveas(gcf, 'voronoi_random.png');
save('vertexPosns_random.mat', 'V');
save('cellsVertices_random.mat', 'C');
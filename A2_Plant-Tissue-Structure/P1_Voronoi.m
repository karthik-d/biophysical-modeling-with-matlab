%% (A) Uniform random points.

N = 50;
randomPosns = rand(2, N);

figure(1);
[V, C, vAll, cAll, xbox, ybox] = getVoronoiDiagram(randomPosns(1, :), randomPosns(2, :), 1);

% save generated data and plot.
saveas(gcf, 'outputs/voronoi_random.fig');
save('outputs/vertexPosns_random.mat', 'V');
save('outputs/cellsVertices_random.mat', 'C');


%% (B)-1 
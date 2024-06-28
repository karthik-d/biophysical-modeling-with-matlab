function PlotConfig(xc, yc, xa, ya, Epara)
%% shrink cancer cell size while increasing cancer cell number
N = Epara.N;
Dc = Epara.Dc;

Na = Epara.Na;
Da = Epara.Da;
D0 = Epara.D0;
idx_start = Epara.idx_start;
idx_end = Epara.idx_end;
xcen = Epara.xcen;
ycen = Epara.ycen;

L = Epara.L;
%%
xshift = [-1; 0; 1; -1; 0; 1; -1; 0; 1] * L;
yshift = [-1; -1; -1; 0; 0; 0; 1; 1; 1] * L;

figure; hold on; box on;

for na = 1:Na
    idx1 = idx_start(na);
    idx2 = idx_end(na);
    for i = 1:9
        if xcen(na) < -0.1 * L - Da(na) || xcen(na) > 1.1 * L + Da(na) ||...
                ycen(na) < -0.1 * L - Da(na) || ycen(na) > 1.1 * L + Da(na)
            continue;
        end
        for ns = idx1:idx2
            rectangle('Position', [xa(ns) + xshift(i) - D0(na) / 2, ya(ns) + yshift(i) - D0(na) / 2, D0(na), D0(na)],...
                      'Curvature', [1 1], 'edgecolor', [150 150 150] / 255, 'facecolor', [150 150 150] / 255);
        end
        patch(xa(idx1:idx2) + xshift(i), ya(idx1:idx2) + yshift(i),...
              [150 150 150] / 255, 'edgecolor', [150 150 150] / 255);
    end
end

if N > 0.5
    x_c_plot = xc;
    y_c_plot = yc;
    x_c_plot = x_c_plot - floor(x_c_plot / L) * L;
    y_c_plot = y_c_plot - floor(y_c_plot / L) * L;
    for i = 1:N
        rad = Dc(i) / 2;
        diam = 2 * rad;
        for j = 1:9
            if x_c_plot(i) < -0.1 * L - rad || x_c_plot(i) > 1.1 * L + rad ||...
                    y_c_plot(i) < -0.1 * L - rad || y_c_plot(i) > 1.1 * L + rad
                continue;
            end
            rectangle('Position',[x_c_plot(i) - rad + xshift(j), y_c_plot(i) - rad + yshift(j), diam, diam],...
                      'Curvature', [1 1], 'FaceColor', [1 0 1], 'Edgecolor', 'k');
        end
    end
end

plot([0 L], [0 0], 'k--', 'linewidth', 1.5)
plot([0 L], [L L], 'k--', 'linewidth', 1.5)
plot([0 0], [0 L], 'k--', 'linewidth', 1.5)
plot([L L], [0 L], 'k--', 'linewidth', 1.5)
axis off
axis equal
xlim([-0.2 1.2] * L)
ylim([-0.2 1.2] * L)
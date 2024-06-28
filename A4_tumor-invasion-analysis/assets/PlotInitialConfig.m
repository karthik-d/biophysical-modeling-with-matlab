function PlotInitialConfig(N)
%% shrink cancer cell size while increasing cancer cell number
amp = 1 / 0.2;
Ncb = floor(N / 2);
Ncs = N - Ncb;
Dc = [1.4 * ones(Ncb, 1); ones(Ncs, 1)];
rsc = 0.5;

KAC = 1; KCC = 1;

Na = 4; Ns_s = 20;
L = sqrt(Na) * amp;
Nb = Na / 2;
Ns_b = round(Ns_s * 1.4);
Ns = [repmat(Ns_b, [Nb, 1]); repmat(Ns_s, [Nb, 1])];
Ns_tot = Nb * (Ns(1) + Ns(Nb+1));

ift = zeros(Ns_tot, 1);
jft = zeros(Ns_tot, 1);
idx_start = zeros(Na,1);
idx_end = zeros(Na,1);

for na = 1:Na
    idx1 = sum(Ns(1:na-1)) + 1;
    idx2 = sum(Ns(1:na));
    idx_start(na) = idx1;
    idx_end(na) = idx2;
    ift(idx1:idx2) = [2:Ns(na) 1]' + idx1 - 1;
    jft(idx1:idx2) = [Ns(na) 1:Ns(na)-1]' + idx1 - 1;
end

a_pos = load('Pos_00001.txt');
Da = amp * a_pos(1) * [1.4 * ones(Nb, 1); ones(Nb, 1)];
D0_ori = Da .* sin(pi ./ Ns);
xy_a_equ = amp * a_pos(2:end);

xa = xy_a_equ(1:Ns_tot);
ya = xy_a_equ(Ns_tot+1:2*Ns_tot);
xa_rsc = zeros(Ns_tot, 1);
ya_rsc = zeros(Ns_tot, 1);

xcen = zeros(Na, 1);
ycen = zeros(Na, 1);

for na = 1:Na
    idx1 = idx_start(na);
    idx2 = idx_end(na);
    xcen(na) = mean(xa(idx1:idx2));
    ycen(na) = mean(ya(idx1:idx2));
    xshift = floor(xcen(na) / L) * L;
    yshift = floor(ycen(na) / L) * L;
    dx = xa(idx1:idx2) - xcen(na);
    dy = ya(idx1:idx2) - ycen(na);
    xa_rsc(idx1:idx2) = rsc * dx + xcen(na) - xshift;
    ya_rsc(idx1:idx2) = rsc * dy + ycen(na) - yshift;
end

xa = xa_rsc;
ya = ya_rsc;
D0 = D0_ori * 1.5;

lx = xa(ift) - xa;
ly = ya(ift) - ya;
lk = sqrt(lx.^2 + ly.^2);

xa_l = zeros(Na, 1);
xa_r = zeros(Na, 1);
ya_u = zeros(Na, 1);
ya_d = zeros(Na, 1);
for na = 1:Na
    idx1 = idx_start(na);
    idx2 = idx_end(na);
    xa_l(na) = min(xa(idx1:idx2));
    xa_r(na) = max(xa(idx1:idx2));
    ya_d(na) = min(ya(idx1:idx2));
    ya_u(na) = max(ya(idx1:idx2));
end

% A = sum(xa .* ya(ift,:) - ya .* xa(ift,:)) / 2;
% lk = sqrt((xa - xa(ift)).^2 + (ya - ya(ift)).^2);
% P0 = zeros(Na, 1);
% for na = 1:Na
%     P0(na) = sum(lk(idx_start(na):idx_end(na)));
% end
% 
% Aa = A + sum(D0 .* P0 / 2 + pi * D0.^2 / 4);
% Aa = A + sum(D0 .* P0 / 2);

AC_para = struct('N', N, 'Na', Na, 'Ns', Ns, 'Ns_tot', Ns_tot,... 
                 'L', L, 'Dc', Dc, 'D0', D0, 'lx', lx, 'ly', ly, 'lk', lk,...
                 'ift', ift, 'jft', jft, 'idx_start', idx_start, 'idx_end', idx_end,...
                 'xa_l', xa_l, 'xa_r', xa_r, 'ya_u', ya_u, 'ya_d', ya_d,...
                 'KAC', KAC, 'KCC', KCC);

Aa = AreaAdipocyte(xa, ya, AC_para);

phi = pi * sum(Dc.^2) / 4 / (L^2 - Aa);
%% randomly generate probe particle positions with particle centers in the void
if N > 0.5
    rng(1);

    xc = zeros(N, 1);
    yc = zeros(N, 1);
    for n = 1:N
        while true
            xc_new = rand(1) * L;
            yc_new = rand(1) * L;
            ol = Overlap(n - 1, xc_new, yc_new, xc, yc, xa, ya, AC_para);
            if ol == 0
                break;
            end
        end
        xc(n) = xc_new;
        yc(n) = yc_new;
    end
    %%
    Fthresh = 1e-10;
    dt_fire = 0.01 * Dc(end);
    Nt_fire = 1e7;
    [xc, yc] = FIRE_VL(xc, yc, xa, ya, AC_para, Fthresh, dt_fire, Nt_fire);
end
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

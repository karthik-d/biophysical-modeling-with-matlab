function ol = Overlap(N_old, xc_new, yc_new, xc, yc, xa, ya, Epara)
%%
ol = 0;

Dc = Epara.Dc;

Na = Epara.Na;
Ns = Epara.Ns;

L = Epara.L;
D0 = Epara.D0;

idx_start = Epara.idx_start;
idx_end = Epara.idx_end;
%%
x_cen = zeros(Na, 1);
y_cen = zeros(Na, 1);
for na = 1:Na
    xcn = xa(idx_start(na):idx_end(na));
    ycn = ya(idx_start(na):idx_end(na));
    x_cen(na) = mean(xcn);
    y_cen(na) = mean(ycn);
end

% particle-obstacle force
for na = 1:Na
    Dnm = (Dc(N_old + 1) + D0(na)) / 2;
    x1 = xc_new - round((xc_new - x_cen(na)) / L) * L;
    y1 = yc_new - round((yc_new - y_cen(na)) / L) * L;
    idx1 = idx_start(na);
    idx2 = idx_end(na);
    xaa = xa(idx1:idx2);
    yaa = ya(idx1:idx2);
    inPoly = InPolyDisk(x1, y1, xaa, yaa, Ns(na));
    if inPoly == 1 % inside of the polygon
        ol = 1;
        return;
    end
    for ns = idx1:idx2
        dx13 = xa(ns) - x1;
        dy13 = ya(ns) - y1;
        d1 = sqrt(dx13^2 + dy13^2);
        if d1 < 0.85 * Dnm
            ol = 1;
            return;
        end
    end
end

% particle-particle force
% for nc = 1:N_old
%     dx = xc_new - xc(nc);
%     dx = dx - round(dx / L) * L;
%     if abs(dx) < 0.85 * Dc
%         dy = yc_new - yc(nc);
%         dy = dy - round(dy / L) * L;
%         if abs(dy) < 0.85 * Dc
%             d = sqrt(dx^2 + dy^2);
%             if d < 0.85 * Dc
%                 ol = 1;
%                 return;
%             end
%         end
%     end
% end

end


function inPoly = InPolyDisk(xp, yp, x_r, y_r, Ns)
%%
inPoly = 1;

ift = [2:Ns 1]';
La = y_r(ift) - y_r;
Lb = x_r - x_r(ift);
Lc = x_r(ift) .* y_r - x_r .* y_r(ift);

xmin = zeros(Ns, 1);
xmax = xmin;
ymin = xmin;
ymax = xmin;
for ns = 1:Ns
    xmin(ns) = min(x_r(ns), x_r(ift(ns)));
    xmax(ns) = max(x_r(ns), x_r(ift(ns)));
    ymin(ns) = min(y_r(ns), y_r(ift(ns)));
    ymax(ns) = max(y_r(ns), y_r(ift(ns)));
end

count = 0;
for ns = 1:Ns
    if abs(La(ns)) < 1e-12
        if xp < xmax(ns)
            count = count + 1;
        end
    else
        if (yp - ymin(ns)) * (yp - ymax(ns)) <= 0
            x0 = -(Lb(ns) * yp + Lc(ns)) / La(ns);
            if x0 > xp
                count = count + 1;
            end
        end
    end
end
if mod(count, 2) == 0
    inPoly = 0;
end

end
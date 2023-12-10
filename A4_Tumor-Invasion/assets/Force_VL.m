function [Fx, Fy] = Force_VL(xc, yc, xa, ya, Epara, VL_c, VL_c_counter, VL_a, VL_a_counter)
%%
Nc = Epara.N;
Dc = Epara.Dc;
KCC = Epara.KCC;

KAC = Epara.KAC;
% KW = Epara.KW;

L = Epara.L;
D0 = Epara.D0;

Fx = zeros(Nc, 1);
Fy = zeros(Nc, 1);
%%
% particle-obstacle force
for vidx = 1:VL_a_counter
    nc = VL_a(vidx, 1);
    na = VL_a(vidx, 2);
    ns = VL_a(vidx, 3);
    Dnm = (Dc(nc) + D0(na)) / 2;
    Dnm_sq = Dnm^2;
    dx13 = xa(ns) - xc(nc);
    dx13 = dx13 - round(dx13 / L) * L;
    if abs(dx13) < Dnm
        dy13 = ya(ns) - yc(nc);
        dy13 = dy13 - round(dy13 / L) * L;
        if abs(dy13) < Dnm
            d1 = sqrt(dx13^2 + dy13^2);
            if d1 < Dnm
                F = KAC * (1 - Dnm / d1) / Dnm_sq;
                Fx(nc) = Fx(nc) + F * dx13;
                Fy(nc) = Fy(nc) + F * dy13;
            end
        end
    end
end

% particle-particle force
for vidx = 1:VL_c_counter
    nc = VL_c(vidx, 1);
    mc = VL_c(vidx, 2);
    Dnm = (Dc(nc) + Dc(mc)) / 2;
    Dnm_sq = Dnm^2;
    dx = xc(mc) - xc(nc);
    dx = dx - round(dx / L) * L;
    if abs(dx) < Dnm
        dy = yc(mc) - yc(nc);
        dy = dy - round(dy / L) * L;
        if abs(dy) < Dnm
            d = sqrt(dx^2 + dy^2);
            if d < Dnm
                F = KCC * (1 - Dnm / d) / Dnm_sq;
                dFx = F * dx;
                dFy = F * dy;
                Fx(nc) = Fx(nc) + dFx;
                Fy(nc) = Fy(nc) + dFy;
                Fx(mc) = Fx(mc) - dFx;
                Fy(mc) = Fy(mc) - dFy;
            end
        end
    end
end
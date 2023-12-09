function [VL_c, VL_c_counter, VL_a, VL_a_counter, xc_save, yc_save] = VerletList(xc, yc, xa, ya, Epara, first_call, VL_c_old, VL_c_counter_old, VL_a_old, VL_a_counter_old, xc_save_old, yc_save_old)
%%
Nc = Epara.N;
Dc = Epara.Dc;

Na = Epara.Na;

L = Epara.L;
D0 = Epara.D0;

idx_start = Epara.idx_start;
idx_end = Epara.idx_end;

r_factor = 1.2;
r_cut = max([Dc; D0(1)]);
r_list = r_factor * r_cut;
r_list_sq = (r_list)^2;
r_skin_sq = ((r_factor - 1) * r_cut)^2;
%%
if first_call < 0.5
    dx = xc - xc_save_old;
    dy = yc - yc_save_old;
    dx = dx - round(dx / L) * L;
    dy = dy - round(dy / L) * L;
    dr = dx.^2 + dy.^2;
    if 4 * max(dr) < r_skin_sq
        VL_c = VL_c_old;
        VL_c_counter = VL_c_counter_old;
        VL_a = VL_a_old;
        VL_a_counter = VL_a_counter_old;
        xc_save = xc_save_old;
        yc_save = yc_save_old;
        return;
    end
end

VL_c = zeros(Nc * 5, 2);
VL_c_counter = 0;

for n = 1:Nc-1
    for m = n+1:Nc
        dx = xc(m) - xc(n);
        dx = dx - round(dx / L) * L;
        if abs(dx) < r_list
            dy = yc(m) - yc(n);
            dy = dy - round(dy / L) * L;
            if abs(dy) < r_list
                dr_sq = dx^2 + dy^2;
                if dr_sq < r_list_sq
                    VL_c_counter = VL_c_counter + 1;
                    VL_c(VL_c_counter, :) = [n, m];
                end
            end
        end
    end
end

VL_a = zeros(Nc * 5, 3);
VL_a_counter = 0;

% particle-obstacle force
for nc = 1:Nc
    x1 = xc(nc);
    y1 = yc(nc);
    for na = 1:Na
        for ns = idx_start(na):idx_end(na)
            dx13 = xa(ns) - x1;
            dx13 = dx13 - round(dx13 / L) * L;
            if abs(dx13) < r_list
                dy13 = ya(ns) - y1;
                dy13 = dy13 - round(dy13 / L) * L;
                if abs(dy13) < r_list
                    d_sq = dx13^2 + dy13^2;
                    if d_sq < r_list_sq
                        VL_a_counter = VL_a_counter + 1;
                        VL_a(VL_a_counter, :) = [nc, na, ns];
                    end
                end
            end
        end
    end
end

xc_save = xc;
yc_save = yc;
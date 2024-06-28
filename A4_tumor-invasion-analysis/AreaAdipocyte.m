function Aa = AreaAdipocyte(xa, ya, Epara)
%%
Na = Epara.Na;
Ns = Epara.Ns;
ift = Epara.ift;
jft = Epara.jft;
idx_start = Epara.idx_start;
idx_end = Epara.idx_end;
lk = Epara.lk;

D0 = Epara.D0;
R0 = D0 / 2;

Ntot = sum(Ns);
%%
xa_inter = zeros(Ntot, 1);
ya_inter = zeros(Ntot, 1);

for na = 1:Na
    R = R0(na);
    R_sq = R^2;
    for ns = idx_start(na):idx_end(na)
        ms = ift(ns);
        d = lk(ns);
        l = lk(ns) / 2;
        h = sqrt(R_sq - l^2);
        xa_inter(ns) = 0.5 * (xa(ms) - xa(ns)) + h / d * (ya(ms) - ya(ns)) + xa(ns);
        ya_inter(ns) = 0.5 * (ya(ms) - ya(ns)) - h / d * (xa(ms) - xa(ns)) + ya(ns);
    end
end

%%
lx_inter = xa_inter - xa_inter(jft);
ly_inter = ya_inter - ya_inter(jft);
lk_inter = sqrt(lx_inter.^2 + ly_inter.^2);

Aa = 0.5 * sum(xa_inter .* ya_inter(ift) - xa_inter(ift) .* ya_inter);

for na = 1:Na
    R = R0(na);
    R_sq = R^2;
    for ns = idx_start(na):idx_end(na)
        theta = asin(lk_inter(ns) / 2 / R);
        Aa = Aa + theta * R_sq - 0.5 * lk_inter(ns) * sqrt(R_sq - (lk_inter(ns) / 2)^2);
    end
end
function [phi, x_t, y_t, vx_t, vy_t, xa, ya, AC_para] = LorentzCancer(N, T)
%% shrink adipocyte and rescale so that the diameter of small cancer cell is 1
amp = 1 / 0.2; % box size rescale factor
Ncb = floor(N / 2);
Ncs = N - Ncb;
Dc = [1.4 * ones(Ncb, 1); ones(Ncs, 1)];
rsc = 0.5; % adipocyte rescale factor

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

a_pos = load('Pos_00001.txt'); % loading adipocyte positions
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

AC_para = struct('N', N, 'Na', Na, 'Ns', Ns, 'Ns_tot', Ns_tot,... 
                 'L', L, 'Dc', Dc, 'D0', D0, 'Da', Da, 'lx', lx, 'ly', ly, 'lk', lk,...
                 'ift', ift, 'jft', jft, 'idx_start', idx_start, 'idx_end', idx_end, 'xcen', xcen, 'ycen', ycen,...
                 'KAC', KAC, 'KCC', KCC);

Aa = AreaAdipocyte(xa, ya, AC_para);

phi = pi * sum(Dc.^2) / 4 / (L^2 - Aa); % packing fraction of cancer cells in the void space
%% randomly generate cancer cell positions with particle centers in the void
rng(1);

xc = zeros(N, 1); % x position for cancer cell center
yc = zeros(N, 1); % y position for cancer cell center
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
%% energy minimize the cancer cell packing
Fthresh = 1e-14;
dt_fire = 0.01;
Nt_fire = 1e7;
[xc, yc] = FIRE_VL(xc, yc, xa, ya, AC_para, Fthresh, dt_fire, Nt_fire);
%% NVT
rng(1);
Nt = 4e4;
dt = 0.01;
beta = 5;

% TO DO
% initialize velocity (variable names: vx and vy) based on the input temperature T here
% vx and vy will both be a N x 1 array
% don't forget to remove the mean velocity in the x and y direction after drawing 
% random numbers to avoid the whole system to drift

vx = normrnd(0, sqrt(T), N, 1);
vy = normrnd(0, sqrt(T), N, 1);

% Remove drift velocity
vx = vx - mean(vx);
vy = vy - mean(vy);

% TO DO
% initialize two N x (Nt/10) arrays, x_t, y_t, to store xc and yc, and two N x (Nt/10)
% arrays, vx_t, vy_t, to store vx and vy, every 10 time steps

x_t = zeros(N, Nt/10);
y_t = zeros(N, Nt/10);

vx_t = zeros(N, Nt/10);
vy_t = zeros(N, Nt/10);

% verlet list
VL_c = zeros(N * 5, 2);
VL_c_counter = 0;
VL_a = zeros(N * 5, 2);
VL_a_counter = 0;
xc_save = xc;
yc_save = yc;
% construct the Verlet list to seped up MD
[VL_c, VL_c_counter, VL_a, VL_a_counter, xc_save, yc_save] = VerletList(xc, yc, xa, ya, AC_para, 1, VL_c, VL_c_counter, VL_a, VL_a_counter, xc_save, yc_save);
% force calculation using the Verlet list
% Fx and Fy are the total x-force and y-force on the cancer cells from
% cancer-cancer and cancer-adipocyte interactions
[Fx, Fy] = Force_VL(xc, yc, xa, ya, AC_para, VL_c, VL_c_counter, VL_a, VL_a_counter);
mass = 1;
k_b = 1;

for nt = 1:Nt
    % TO DO
    % update velocity and position using the "BAOAB" method here
    % update velocity for half step
	vx_half = vx + (dt/2)*Fx/mass;
	vy_half = vy + (dt/2)*Fy/mass;    
    % update position for half step
    rx_half = xc + (dt/2)*vx;
	ry_half = yc + (dt/2)*vy;
    % update velocity due to thermal noise
	vx_rate_half = exp(-1*(beta*dt)/mass)*vx_half + sqrt(1-exp(-2*(beta*dt)/mass))*sqrt(k_b*T/mass)*randn;
	vy_rate_half = exp(-1*(beta*dt)/mass)*vy_half + sqrt(1-exp(-2*(beta*dt)/mass))*sqrt(k_b*T/mass)*randn;
    % update position for half step
	xc = rx_half + (dt/2)*vx_rate_half;
	yc = ry_half + (dt/2)*vy_rate_half;
	if mod(nt, 10) == 0
		x_t(:, nt/10) = xc;
		y_t(:, nt/10) = yc;
	end
    % calculate updated force after a full step
    [VL_c, VL_c_counter, VL_a, VL_a_counter, xc_save, yc_save] = VerletList(xc, yc, xa, ya, AC_para, 0, VL_c, VL_c_counter, VL_a, VL_a_counter, xc_save, yc_save);
    [Fx, Fy] = Force_VL(xc, yc, xa, ya, AC_para, VL_c, VL_c_counter, VL_a, VL_a_counter);
    % update velocity for half step
    vx = vx_rate_half + ((dt/2)*Fx)/mass;     % still, Fi=0.
	vy = vy_rate_half + ((dt/2)*Fy)/mass;
	if mod(nt, 10) == 0
		vx_t(:, nt/10) = vx;
		vy_t(:, nt/10) = vy;
	end
end
%%

function [x_equ, y_equ] = FIRE_VL(x, y, xa, ya, Epara, Fthresh, dt_fire, Nt)
%% This is based on FIRE 2.0 provided in https://arxiv.org/pdf/1908.02038.pdf
%% Velocity verlet is chosen as the MD simulator
%%
N = Epara.N;
%% FIRE parameter for fast energy minimization
N_delay = 20;
N_pn_max = 2000;
f_inc = 1.1;
f_dec = 0.5;
a_start = 0.15;
f_a = 0.99;
dt_max = 10 * dt_fire;
dt_min = 0.05 * dt_fire;
initialdelay = 1;
%% Recording
% F2norm = zeros(Nt, 1);
% Eval = zeros(Nt, 1);
%% VL initialization
VL_c = zeros(N * 5, 2);
VL_c_counter = 0;
VL_a = zeros(N * 5, 2);
VL_a_counter = 0;
xc_save = x;
yc_save = y;
[VL_c, VL_c_counter, VL_a, VL_a_counter, xc_save, yc_save] = VerletList(x, y, xa, ya, Epara, 1, VL_c, VL_c_counter, VL_a, VL_a_counter, xc_save, yc_save);
%% Initialization
a_fire = a_start;
delta_a_fire = 1 - a_fire;
dt = dt_fire;
dt_half = dt / 2;
Vx = zeros(N, 1);
Vy = zeros(N, 1);
[Ax, Ay] = Force_VL(x, y, xa, ya, Epara, VL_c, VL_c_counter, VL_a, VL_a_counter);
if max(abs([Ax; Ay])) < Fthresh
    x_equ = x;
    y_equ = y;
    return;
end

N_pp = 0; % number of P being positive
N_pn = 0; % number of P being negative
it = 0;
%% FIRE
while (true)
    it = it + 1;
    % FIRE update
    P = Ax' * Vx + Ay' * Vy;
    
    if P > 0
        N_pp = N_pp + 1;
        N_pn = 0;
        if N_pp > N_delay
            dt = min(f_inc * dt, dt_max);
            dt_half = dt / 2;
            a_fire = f_a * a_fire;
            delta_a_fire = 1 - a_fire;
        end
    else
        N_pp = 0;
        N_pn = N_pn + 1;
        if N_pn > N_pn_max
            break;
        end
        if (initialdelay < 0.5) || (it >= N_delay)
            if f_dec * dt > dt_min
                dt = f_dec * dt;
                dt_half = dt / 2;
            end
            a_fire = a_start;
            delta_a_fire = 1 - a_fire;
            x = x - Vx * dt_half;
            y = y - Vy * dt_half;
            Vx = zeros(N, 1);
            Vy = zeros(N, 1);
        end
    end
    
    % MD using Verlet method
    Vx = Vx + Ax * dt_half;
    Vy = Vy + Ay * dt_half;
    rsc = a_fire * norm([Vx; Vy]) / norm([Ax; Ay]);
    Vx = delta_a_fire * Vx + rsc * Ax;
    Vy = delta_a_fire * Vy + rsc * Ay;
    x = x + Vx * dt;
    y = y + Vy * dt;
    [VL_c, VL_c_counter, VL_a, VL_a_counter, xc_save, yc_save] = VerletList(x, y, xa, ya, Epara, 0, VL_c, VL_c_counter, VL_a, VL_a_counter, xc_save, yc_save);
    [Ax, Ay] = Force_VL(x, y, xa, ya, Epara, VL_c, VL_c_counter, VL_a, VL_a_counter);
    Vx = Vx + Ax * dt_half;
    Vy = Vy + Ay * dt_half;
    
    % For recording
%     Eval(it) = E;
%     F2norm(it) = sum(grad.^2);

    % exit when forces are smaller than threshold
    if max(abs([Ax; Ay])) < Fthresh || it > Nt
        break;
    end
end
%%
x_equ = x;
y_equ = y;
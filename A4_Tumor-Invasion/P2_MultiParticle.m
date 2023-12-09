
N = [20 42 43 44];
T = 10e-12;

N = 44;
[phi, x_t, y_t, vx_t, vy_t, xa, ya, AC_para] = LorentzCancer(N, T);

T_inst = (vx_t.^2 + vy_t.^2)/(2*N);

figure(1);
plot(length(T_inst), T_inst);
hold off;
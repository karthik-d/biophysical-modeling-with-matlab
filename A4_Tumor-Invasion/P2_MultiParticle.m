N_v = [20 42 43 44];
T = 10e-12;

for i = 1:length(N_v)
	N = N_v(i);
	[phi, x_t, y_t, vx_t, vy_t, xa, ya, AC_para] = LorentzCancer(N, T);
	% plot for N=44.
	if N==44
		T_inst = sum(vx_t.^2 + vy_t.^2)/(2*N);
		disp(mean(T_inst));

		figure(1);
		plot(1:length(T_inst), T_inst);
		hold off;
	end
end


% the logspace trick
x_samples = round(logspace(0, 6, 500), 0);
MSD = zeros(length(x_samples), length(N_particles));
deltaT = zeros(length(x_samples), length(N_particles));
% N = 44
for i = 1:length(x_samples)
    idx = x_samples(i);
    delta_x = x_t(:, 1+idx:end) - x_t(:, 1:end-idx);
    delta_y = y_t(:, 1+idx:end) - y_t(:, 1:end-idx);
    MSD(i,1) = mean(delta_x.^2+delta_y.^2, 'all');
	deltaT(i, 1) = dt*i;
end


% ideal exponential decay.
expected = exp(-5*deltaT);

figure(2);
% plot figures.
plot(deltaT(:, 1), MSD(:, 1), "-");
plot(deltaT(:, 1), expected(:, 1), ":");
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
hold on;
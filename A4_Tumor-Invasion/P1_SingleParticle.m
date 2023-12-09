N_particles = 1;
% PlotInitialConfig(N_particles);

N_timesteps = 1e5;
T_step = 1e-4;

mass = 1;
k_b = 1;

% initialize positions.
rx = zeros(N_particles, N_timesteps);
ry = zeros(N_particles, N_timesteps);

% initialize velocities.
vx = randn(N_particles, N_timesteps);
vy = randn(N_particles, N_timesteps);

% compute temperature - using equipartition theorem.
T_current = sum(vx.^2 + vy.^2)/(2*N_particles);
T_target = 1e-4;  

% scale velocities to achieve target temperature.
scale_factor = sqrt(T_target/T_current);
vx = vx .* scale_factor;
vy = vy .* scale_factor;

% drag values.
beta_drag_v = [2 10 50 100];
Fi = zeros(N_particles, N_timesteps);

% init figure.
figure(1), clf, hold on, box on;

for bi = 1:length(beta_drag_v)

	T = zeros(N_particles, N_timesteps);
	ti = 1;
	for i = 1:(N_timesteps-1)
		t_new = ti + T_step;
		i_new = i + 1;
		% update half-velocity.
		vx_half = vx(1, i) + (T_step/2)*Fi(1, i)/mass;
		vy_half = vy(1, i) + (T_step/2)*Fi(1, i)/mass;
		% update half-posn.
		rx_half = rx(1, i) + (T_step/2)*vx(1, i);
		ry_half = ry(1, i) + (T_step/2)*vy(1, i);
		% update half-velocity rate.
		vx_rate_half = exp(-1*(beta_drag_v(bi)*T_step)/mass)*vx_half + sqrt(1-exp(-2*(beta_drag_v(bi)*T_step)/mass))*sqrt(k_b*T_target/mass)*randn;
		vy_rate_half = exp(-1*(beta_drag_v(bi)*T_step)/mass)*vy_half + sqrt(1-exp(-2*(beta_drag_v(bi)*T_step)/mass))*sqrt(k_b*T_target/mass)*randn;
		% update posn.
		rx(1, i_new) = rx_half + (T_step/2)*vx_rate_half;
		ry(1, i_new) = ry_half + (T_step/2)*vy_rate_half;
		% update velocity.
		vx(1, i_new) = vx_rate_half + ((T_step/2)*Fi(1, i_new))/mass;     % still, Fi=0.
		vy(1, i_new) = vy_rate_half + ((T_step/2)*Fi(1, i_new))/mass;
		% update time.
		ti = t_new;
		T(1, i_new) = t_new;
		% compute autocorrelation.
	end

	% reference plot.
	auto_corrs = zeros(N_timesteps, 1);
	deltaT = zeros(N_timesteps, 1);
	i = 1;
	for t = 1:(N_timesteps-1)
		window = t;
		v_numerator = mean(dot(vx(1, 1+window:end), vx(1, 1:end-window)) + dot(vy(1, 1+window:end), vy(1, 1:end-window)));
		v_denominator = mean(dot(vx(1, 1+window:end), vx(1, 1+window:end)) + dot(vy(1, 1+window:end), vy(1, 1+window:end)));
		auto_corrs(i) = v_numerator / v_denominator;
		deltaT(i) = t*T_step;
		i = i + 1;
	end

	% ideal exponential decay.
	expected = exp(-1*beta_drag_v(bi)*deltaT/mass);

	% plot figures.
	ax = gca;
	ax.XScale = 'log';
	plot(deltaT, auto_corrs, "-");
	plot(deltaT, expected, ":");
	hold on;
end

legend('Actual $\beta=2$', 'Ideal $\beta=2$', ...
       'Actual $\beta=10$', 'Ideal $\beta=10$', ...
       'Actual $\beta=50$', 'Ideal $\beta=50$', ...
       'Actual $\beta=100$', 'Ideal $\beta=100$', ...
       'Interpreter', 'Latex', 'Location','southwest');


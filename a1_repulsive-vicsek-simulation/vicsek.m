function [xtotal, ytotal, pols] = vicsek(N,pf,r0,v0,dt,eta,beta,Nsteps,Nplot)
%% Test function to loop over vicsek simulation, output trajectories

% Box Size
L = 1.0;

% Definte particle sizes (rc) based on packing fraction
rc = 2.0*L*sqrt(pf/(pi*N));

% scale r0 by rc
r0 = r0*rc;

% initialize trajectory values
xtotal = zeros(N,Nsteps);
ytotal = zeros(N,Nsteps);

% initialize polarization values
pols = zeros(Nsteps,1);

% Initialize Positions and Velocities randomly
rs = rand(N, 2) * L;
vs = randn(N, 2);
vnorm = sqrt(sum(vs'.^2))';
vs = vs .* v0 ./ [vnorm, vnorm];

% Loop over simulation
for step = 1:Nsteps
    % Perform Integration
    vs = vicsekvelocities(N,v0,r0,rc,eta,beta,L,rs,vs) ;
    rs = rs + vs * dt;
    
    % store polarization
    meanVs = mean(vs);
    meanVsNrm = sqrt(meanVs(1)^2 + meanVs(2)^2);
    pols(step) = meanVsNrm/v0;
    
    % store trajectories
    xtotal(:,step) = rs(:,1);
    ytotal(:,step) = rs(:,2);
    
    % Plot cells and velocities
    if mod(step, Nplot) == 0
        fprintf('On step %d for beta = %f, eta = %f\n',step,beta,eta);
    end
end



end
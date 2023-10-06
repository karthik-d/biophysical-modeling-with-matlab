function []=vicsekplot(rs, vs, rc, L)
    %{
    Function that plots disks of diameter rc in a box of length L. 
    Plots vectors indicating the direction of the velocities on top of the disks.
    
    Inputs:
        rs         Nx2 matrix of disk positions in the box, where N is the number of disks
        vs         Nx2 matrix of disk velocities
        rc         Diameter of the disks
        L          Length of the side of the square box
        
    Outputs:
        Overwrites the figure plotted with the function handle, f_h
    %}
    
    % Calculate the number of particles
    [N, ~] = size(rs);
    
    % Clear the figure
    clf
    
    % Plot velocities as a quiver plot
    quiver(mod(rs(:,1),L),mod(rs(:,2),L),vs(:,1),vs(:,2),0.25)
    
    % Define the plot edges to be the edges of the box
    xlim([0, L]);
    ylim([0, L]);
    
    % Plot circles with proper radii for the disks
    h=zeros(1,N);
    for ind=1:N
        h(ind)=rectangle('Position',[mod(rs(ind,1)-.5*rc,L) mod(rs(ind,2)-.5*rc,L) rc rc],'Curvature',[1 1],'edgecolor','b');
    end
    drawnow;
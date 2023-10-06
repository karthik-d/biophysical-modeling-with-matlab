function vs = vicsekvelocities(N, v0, r0, rc, eta, beta, L, rs, vs)
    %{
    Function that finds the updated velocity for a Vicsek model with repulsion.
    The model is from PhysRevE.77.046113
    
    Inputs:
        N          The number of cells
        v0         The constant speed of the cells
        r0         The distance cutoff within which cells are considered neighbors 
        rc         The diameter of the cell
        eta        The relative strength of the noise
        beta       The relative strength of the repulsion
        L          The size of the box
        rs         Nx2 matrix of cell positions
        vs         Nx2 matrix of cell velocities
    
    Outputs:
        vs     Nx2 matrix of updated velocities due to the repulsive Vicsek algorithm
    
    %}
    
    % Initialize Variables
    sum_vs = ;   % Sum of neighboring velocities
    Si_norm = ;  % Number of neighbors for each cell i
    Fi = ;       % Total repulsive force for each cell i
    
    % Calculate Fi
    for i=1:N
       % Calculate distance between cell i and the rest within a periodic box
       rijs = -[rs(:, 1) - rs(i, 1), rs(:, 2) - rs(i, 2)];
       rijs = mod((rijs + L./2), L) - L./2;
       dists = sqrt(sum(rijs'.^2))';
       
       % Calculate the set Si 
       Si = ;
       Si_norm(i) = ;
       
       % Calculate sum_vs(i,:)
       sum_vs(i, :) = ;
       
       % Calculate the repulsive force due to each cell on i
       
       Fi(i,:)=;
    end
    
    % Calculate a matrix of random unit vectors
    noise = normer(randn(N, 2));
    
    vs = normer(sum_vs./v0 + beta.*Fi + eta.*[Si_norm,Si_norm].*noise) .* v0;
end


function v = normer(input)
    %{
    A function that divides each row in a matrix by its norm.
    
    Inputs:
        input     A matrix for which to normalize along the second axis
    Outputs:
        v         The normalized matrix
    %}
    
    normval = sqrt(sum(input.^2, 2)) + 10^-16;
    v = input ./ [normval, normval];
end
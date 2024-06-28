function area = triangularArea(p1,p2,p3)
%% FUNCTION to calculate area of triangle given three points in 2D

% check that arguments are the right size
if sum(size(p1) ~= [1,2]) == 2 && sum(size(p1) ~= [2,1]) == 2
    error('Argument for first point in input is the wrong size\n');
end

if sum(size(p2) ~= [1,2]) == 2 && sum(size(p2) ~= [2,1]) == 2
    error('Argument for second point in input is the wrong size\n');
end

if sum(size(p3) ~= [1,2]) == 2 && sum(size(p3) ~= [2,1]) == 2
    error('Argument for third point in input is the wrong size\n');
end

% calculate distance vectors
V1 = p2 - p1;
V2 = p3 - p1;

% calculate area
area = 0.5*abs(V1(1)*V2(2) - V2(1)*V1(2));

end
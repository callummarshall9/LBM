function index = scalar_index(x,y,z,NX,NY,NZ)
%SCALAR_INDEX Summary of this function goes here
%   Detailed explanation goes here
x_index = x - 1;
y_index = y - 1;
z_index = z - 1;
index = (x_index + y_index * NX + z_index * NY)+1;
end


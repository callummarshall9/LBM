function visualise_density(actual_name,NX,NY,NZ,index)
M = csvread(actual_name, 1,0);
size = NX * NY * NZ;
density = M(1:size,1);
density_3d = zeros(NX,NY,NZ);
for i=1:NX
    for j=1:NY
        for k=1:NZ
            density_3d(i,j,k)= density(scalar_index(i,j,k,NX,NY,NZ));
        end
    end
end
x = [0:1:NX-1];
y = [0:1:NY-1];
z = [0:1:NZ-1];

xslice = [0:5:(NX-1)];   
yslice = [];
zslice = [];
density_title = strcat("Density plot of LBM (",actual_name,")");
fig = figure('Name',density_title);
slice(x,y,z,density_3d,xslice,yslice,zslice)
saveas(fig,strcat(actual_name,".png"));
end

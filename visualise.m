fname = 'options.json';
val = jsondecode(fileread(fname));
grid_size = val.grid_size;
NX=grid_size(1);
NY=grid_size(2);
NZ=grid_size(3);
save_every = val.save_every;
indices=[0:save_every:1000];
indice_size = length(indices);
for i = 1:indice_size
    indice = indices(i);
    csv_name = strcat("output/",string(indice),".csv");
    visualise_density(csv_name,NX,NY,NZ,i);
end
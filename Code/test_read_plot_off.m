name = 'Data/shrec10/0003.null.0.off';
[vertex,faces] = read_off(name);
% display the mesh
clf;
plot_mesh(vertex, faces);
shading interp;
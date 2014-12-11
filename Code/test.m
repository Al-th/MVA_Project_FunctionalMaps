%clear; clc;
init;

%%
%load shape 1
name = 'Data/shrec10/0003.null.0.off';
shape1.name = name;
[vertex,faces] = read_off(name);
vertex = vertex'; faces = faces';
shape1.vertex = vertex;
shape1.faces = faces;
%Compute ALB
[PHI,E,L,Am] = ALB_spectrum(vertex,faces);
shape1.phi = PHI;
shape1.E = E;
shape1.Am = Am;
%Compute HKS
hks = HKS(PHI, E, diag(Am), false);
shape1.HKS = hks;
%Compute WKS
WKS = compute_WKS_from_EV(E,PHI);
shape1.WKS = WKS;
clear WKS; clear hks;
clear PHI; clear E; clear L; clear Am;
clear vertex; clear faces; clear name;

%%
%load shape 2
name = 'Data/shrec10/0003.isometry.2.off';
shape2.name = name;
[vertex,faces] = read_off(name);
vertex = vertex'; faces = faces';
shape2.vertex = vertex;
shape2.faces = faces;
[PHI,E,L,Am] = ALB_spectrum(vertex,faces);
shape2.phi = PHI;
shape2.E = E;
shape2.Am = Am;
%Compute HKS
hks = HKS(PHI, E, diag(Am), false);
shape2.HKS = hks;
%Compute WKS
WKS = compute_WKS_from_EV(E,PHI);
shape2.WKS = WKS;
clear WKS; clear hks;
clear PHI; clear E; clear L; clear Am;
clear vertex; clear faces; clear name;

%%
%Compute connected component of shape 1

[C,PD] = compute_connected_component(shape1,80);
figure(1);
options.face_vertex_color = C;
plot_mesh(shape1.vertex,shape1.faces,options);
shading interp; colormap jet(256);
%%
%Compute connected component of shape 2

[C2,PD] = compute_connected_component(shape2,100);
figure(2);
options.face_vertex_color = C2;
plot_mesh(shape2.vertex,shape2.faces,options);
shading interp; colormap jet(256);



%%
% Add constraint to linear system
a1 = shape1.phi'*shape1.Am*shape1.HKS;
a2 = shape1.wks_phi'*shape1.WKS;
a = [a1 a2];

b1 = shape2.phi'*shape2.Am*shape2.HKS;
b2 = shape2.wks_phi'*shape2.WKS;
b = [b1 b2];




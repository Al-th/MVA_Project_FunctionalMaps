% Test compute ALB HKS and WKS

%% Load mesh
name = 'Data/shrec10/0003.isometry.5.off';
[vertex,faces] = read_off(name);
faces = faces'; vertex = vertex';

%% Compute Anisotropic Laplace-beltrami spectrum
[PHI,E,L,Am] = ALB_spectrum(vertex,faces);

%% Compute Wave Kernel Signature
[WKS,E_wks,PHI_wks,L_wks] = compute_wks(vertex,faces);

%% Compute Heat Kernel Signature

hks = HKS(PHI, E, diag(Am), false);

%%
clear; clc;
init;

%load shape 1
name = 'Data/shrec10/0003.null.0.off';
shape1.name = name;
[vertex,faces] = read_off(name);
vertex = vertex'; faces = faces';
shape1.vertex = vertex;
shape1.faces = faces;
[PHI,E,L,Am] = ALB_spectrum(vertex,faces);
shape1.phi = PHI;
shape1.E = E;
hks = HKS(PHI, E, diag(Am), false);
shape1.HKS = hks;
clear hks;
clear PHI; clear E; clear L; clear Am;
clear vertex; clear faces; clear name;

%load shape 2
name = 'Data/shrec10/0003.isometry.5.off';
shape2.name = name;
[vertex,faces] = read_off(name);
vertex = vertex'; faces = faces';
shape2.vertex = vertex;
shape2.faces = faces;
[PHI,E,L,Am] = ALB_spectrum(vertex,faces);
shape2.phi = PHI;
shape2.E = E;
hks = HKS(PHI, E, diag(Am), false);
shape2.HKS = hks;
clear hks;
clear PHI; clear E; clear L; clear Am;
clear vertex; clear faces; clear name;



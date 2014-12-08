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
%Compute ALB
[PHI,E,L,Am] = ALB_spectrum(vertex,faces);
shape1.phi = PHI;
shape1.E = E;
shape1.Am = Am;
%Compute HKS
hks = HKS(PHI, E, diag(Am), false);
shape1.HKS = hks;
%Compute WKS
[WKS,E_wks,PHI_wks,L_wks] = compute_wks(vertex,faces);
shape1.WKS = WKS;
shape1.wks_phi = PHI_wks;
shape1.wks_E = E_wks;
clear WKS; clear E_wks; clear PHI_wks; clear L_wks;
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
shape2.Am = Am;
hks = HKS(PHI, E, diag(Am), false);
shape2.HKS = hks;
[WKS,E_wks,PHI_wks,L_wks] = compute_wks(vertex,faces);
shape2.WKS = WKS;
shape2.wks_phi = PHI_wks;
shape2.wks_E = E_wks;
clear WKS; clear E_wks; clear PHI_wks; clear L_wks;
clear hks;
clear PHI; clear E; clear L; clear Am;
clear vertex; clear faces; clear name;


a1 = shape1.phi(1:10000,:)'*shape1.HKS(1:10000,:);
a2 = shape1.phi(10001:end,:)'*shape1.HKS(10001:end,:);
a3 = shape1.wks_phi'*shape1.WKS;
a = [a1 a2 a3];

b1 = shape2.phi(1:10000,:)'*shape2.Am*shape2.HKS(1:10000,:);
b2 = shape2.wks_phi(10001:end,:)'*shape2.WKS(10001:end,:);
b = [b1 b2];




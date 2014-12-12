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
name = 'Data/shrec10/0003.isometry.1.off';
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

C1 = persistance_based_segmentation(shape1,7);
shape1.connected_component = C1;
list_label_C1 = union(C1,C1);

% figure(1);
% options.face_vertex_color = compute_color_from_connected_component(C1, list_label_C1);
% plot_mesh(shape1.vertex,shape1.faces,options);
% shading interp; colormap jet(256);
%%
%Compute connected component of shape 2

C2 = persistance_based_segmentation(shape2,7);
shape2.connected_component = C2;
list_label_C2 = union(C2,C2);

% figure(2);
% options.face_vertex_color = compute_color_from_connected_component(C2,list_label_C2);
% plot_mesh(shape2.vertex,shape2.faces,options);
% shading interp; colormap jet(256);

%%
%Matching segment bewteen shape1 and shape2
list_descriptors_C1 = compute_descriptors(shape1);

list_descriptors_C2 = compute_descriptors(shape2);

[~,perm1] = sort(list_descriptors_C1);
[~,perm2] = sort(list_descriptors_C2);

list_matching = [];
for i1=1:size(list_descriptors_C1,1)
    %Find best match on list_descriptors_C2
    desc1 = list_descriptors_C1(i1);
    for i2=1:size(list_descriptors_C2,1)
        
    end
end

diff = abs(list_descriptors_C1(perm1) - list_descriptors_C2(perm2))...
    ./(list_descriptors_C1(perm1) + list_descriptors_C2(perm2));

list_matching = [ list_label_C1(perm1) list_label_C2(perm2) diff];

%Plot mesh with color correspondence to same segment
figure(1);
options.face_vertex_color = compute_color_from_connected_component(C1, list_matching(:,1));
plot_mesh(shape1.vertex,shape1.faces,options);
shading interp; colormap jet(256);

figure(2);
options.face_vertex_color = compute_color_from_connected_component(C2, list_matching(:,2));
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




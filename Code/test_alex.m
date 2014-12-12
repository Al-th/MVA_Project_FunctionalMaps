clear; clc;
init;
%%

alboptions.n_eigenvalues = 100; 
%%
%load shape 1
name = 'Data/shrec10/0003.null.0.off';
shape1.name = name;
[vertex,faces] = read_off(name);
vertex = vertex'; faces = faces';
shape1.vertex = vertex;
shape1.faces = faces;
%Compute ALB
[PHI,E,L,Am] = ALB_spectrum(vertex,faces,alboptions);
shape1.phi = PHI;
shape1.E = E;
shape1.Am = Am;
shape1.L  = L;
%Compute HKS
hks = HKS(PHI, E, diag(Am), false);
shape1.HKS = hks;
%Compute WKS
WKS = compute_WKS_from_EV(E,PHI,alboptions);
shape1.WKS = WKS;
clear WKS; clear hks;
clear PHI; clear E; clear L; clear Am;
clear vertex; clear faces; clear name;

%%
%load shape 2
name = 'Data/shrec10/0003.isometry.5.off';
shape2.name = name;
[vertex,faces] = read_off(name);
vertex = vertex'; faces = faces';
shape2.vertex = vertex;
shape2.faces = faces;
[PHI,E,L,Am] = ALB_spectrum(vertex,faces,alboptions);
shape2.phi = PHI;
shape2.E = E;
shape2.Am = Am;
shape2.L = L;
%Compute HKS
hks = HKS(PHI, E, diag(Am), false);
shape2.HKS = hks;
%Compute WKS
WKS = compute_WKS_from_EV(E,PHI,alboptions);
shape2.WKS = WKS;
clear WKS; clear hks;
clear PHI; clear E; clear L; clear Am;
clear vertex; clear faces; clear name;


%%
%Test KD-tree
tree1 = kdtree_build(shape1.phi);
tree2 = kdtree_build(shape2.phi);

%
clear options
colors = zeros(19248,1);
colors(1:1000,1) = 2;
options.face_vertex_color = colors
plot_mesh(shape1.vertex,shape1.faces,options);
shading interp; colormap jet(256);

%%
%Search for 100 first vertex and set colors 
clear options2;
color = zeros(19248,1);
for i = 1:1000
    p1 = shape1.phi(i,:)';
    nn = kdtree_k_nearest_neighbors(tree2,C*p1,1);
    [fi,fj] = find(shape2.faces==nn);
    color(shape2.faces(fi,1)) = 2;
    color(shape2.faces(fi,2)) = 2;
    color(shape2.faces(fi,3)) = 2;
end
options2.face_vertex_color = color;
plot_mesh(shape2.vertex,shape2.faces,options2);
shading interp; colormap jet(256);


%%
%How many eigenfunctions are needed
%Based on the idea of
%http://alice.loria.fr/publications/papers/2006/SMI_Laplacian/SMI_Laplacian.pdf
%to reconstruct functions based on projections
%
% 1. Project function on the eigenfunctions of the LB operator
% 2. Reconstruct function based on first eigenvalues
% 3. Take the sum of the function to see how it converges
func = shape1.WKS(:,5)
shape1.projectedFunction = shape1.phi'*shape1.Am*func;
shape1.reconstructedFunction = zeros(19248,1);

s = [];
for i = 1:100
    shape1.reconstructedFunction = shape1.reconstructedFunction + shape1.projectedFunction(i)*shape1.phi(:,i);
    s = [s sum(sum(shape1.reconstructedFunction))/19248];
    figure(1);
    clf
    options.face_vertex_color = shape1.reconstructedFunction;
    plot_mesh(shape1.vertex,shape1.faces,options);
    shading interp; colormap jet(256);
    pause(0.1);
end
plot(s);
hold on
plot(repmat(sum(sum(func))/19248,1,100),'r')
%%
figure(1);
options.face_vertex_color = shape1.HKS(:,2);
plot_mesh(shape1.vertex,shape1.faces,options);
shading interp; colormap jet(256);

%%
% Post processing iterative refinement
tree2 = kdtree_build(shape2.phi);
C = eye(100,100);

C_iteration = C;

for iteration = 1:1000
    CPHI_M = C_iteration*shape1.phi';
    correspondences = zeros(19248,1);
    for i = 1:19248
        if(mod(i,100)==0)
            i
        end
        [idx] = kdtree_k_nearest_neighbors(tree2,CPHI_M(:,i),1);
        correspondences(i) = idx;
    end
    e = CPHI_M - shape2.phi(correspondences,:);
    sum(sum(e))
end

%% Writing linear constraints for C determination

A = [];
b = [];

%%add descriptor preservation constraints
for i = 1:size(shape1.HKS,2)
    [dA,db] = computeDescriptorConstraints(shape1.HKS(:,i),shape1,shape2.HKS(:,i),shape2);
    A = sparse([A ; dA]);
    b = sparse([b ; db]);
end

%%
for i = 1:size(shape1.WKS,2)
    [dA,db] = computeDescriptorConstraints(shape1.WKS(:,i),shape1,shape2.WKS(:,i),shape2);
    A = sparse([A ; dA]);
    b = sparse([b ; db]);
end

%%add operator commutativity constraints
for i = 1:size(shape1.HKS,2)
    [dA,db] = computeOperatorCommutativityConstraints(shape1.HKS(:,i),shape1,shape2);
    A = sparse([A; dA]);
    b = sparse([b; db]);
end

for i = 1:size(shape1.WKS,2)
    [dA,db] = computeOperatorCommutativityConstraints(shape1.WKS(:,i),shape1,shape2);
    A = sparse([A; dA]);
    b = sparse([b; db]);
end



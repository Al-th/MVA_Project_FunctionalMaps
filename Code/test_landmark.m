clear; 
clc;
init;

%% Loading the datas
name1 = 'Data/shrec10/0003.null.0.off';
name2 = 'Data/shrec10/0003.isometry.5.off';
gt = load('Data\shrec10gt\0003.isometry.5.labels');
%%Avec la 4 c'est top horse
%%Avec la 5 c'est top horse
%%Avec la 2 c'est top dog

shape1 = getShape(name1);
shape2 = getShape(name2);

%% Find connected components

C1 = persistance_based_segmentation(shape1,7);
shape1.connected_component = C1;
list_label_C1 = union(C1,C1);


C2 = persistance_based_segmentation(shape2,7);
shape2.connected_component = C2;
list_label_C2 = union(C2,C2);

list_matching = [];

list_descriptors_C1 = compute_descriptors_for_matching(shape1);
list_descriptors_C2 = compute_descriptors_for_matching(shape2);

nb_comp_C1 = size(list_label_C1,1);
nb_comp_C2 = size(list_label_C2,1);

err=[];
for i=1:nb_comp_C1
    for j=1:nb_comp_C2
        err(i,j) = (list_descriptors_C1(i) - list_descriptors_C2(j))^2./(list_descriptors_C1(i) + list_descriptors_C2(j));
    end
end
i = 0;
INFTY = max(max(err))+1;
while(i~=min(nb_comp_C1,nb_comp_C2))
    i = i+1;
    [min_per_col,idx_row] = min(err);
    [diff,idx_col] = min(min_per_col);
    list_matching(i,:) = [ list_label_C1(idx_row(idx_col)) list_label_C2(idx_col) diff ];
    
    % Replace idx_row(idx_col)th row and idx_colth col of err with a high
    % value
    err(idx_row(idx_col),:) = INFTY;
    err(:,idx_col) = INFTY;
end

%% Find index corresponding to head an tail on each shape

[test, index] = sort(shape1.HKS(C1==11369,1),'descend');
%[my, idy] = max(shape1.vertex(:,1));
[M1x, s1_tail] = max(shape1.vertex(:,2));
[m1x, s1_head] = min(shape1.vertex(:,2));
[M2x, s2_tail] = max(shape2.vertex(:,2));
[m, s2_head] = max(1.*(gt==11722));

%% Check landmark position
figure();

options.face_vertex_color = compute_color_from_connected_component(C2,10392 );
%options.view_param = [180,90];
plot_mesh(shape2.vertex,shape2.faces,options);
hold on;
plot3(shape2.vertex(s2_head,1),shape2.vertex(s2_head,2),shape2.vertex(s2_head,3),'or')
plot3(shape2.vertex(s2_tail,1),shape2.vertex(s2_tail,2),shape2.vertex(s2_tail,3),'or')
shading interp; colormap jet(256);
hold off;

%% Create parts indicator functions

shape1.parts = [];
shape2.parts = [];

% shape1.parts = [shape1.parts 1.*(C1==15346)];
% shape2.parts = [shape2.parts getAssociatedSegmentFromNull(1.*(C1==15346),gt)];
% shape1.parts = [shape1.parts 1.*(C1==11369)];
% shape2.parts = [shape2.parts getAssociatedSegmentFromNull(1.*(C1==11369),gt)];
% shape1.parts = [shape1.parts 1.*(C1==4833)];
% shape2.parts = [shape2.parts getAssociatedSegmentFromNull(1.*(C1==4833),gt)];

%% Add constraint for landmark (exactly like segment)
tmp = 1:numel(shape1.vertex(:,1));
tmp = tmp';
[D1,S,Q] = perform_fast_marching_mesh(shape1.vertex,shape1.faces,s1_tail);
[D2,S,Q] = perform_fast_marching_mesh(shape2.vertex,shape2.faces,s2_head);
[D3,S,Q] = perform_fast_marching_mesh(shape1.vertex,shape1.faces,s1_head);
[D4,S,Q] = perform_fast_marching_mesh(shape2.vertex,shape2.faces,s2_tail);

shape1.parts = [shape1.parts D1];
shape2.parts = [shape2.parts D2];
shape1.parts = [shape1.parts D3];
shape2.parts = [shape2.parts D4];

%% Create constraints
shape1.fun_segment = [];
shape2.fun_segment = [];

W = zeros(1206,1);
weightWKS = 0.01;
weightHKS = 0.01;
weightSegment = 1;
weightWKSComm = 0.01;
weightHKSComm = 0.01;
weightSegmentComm = 1;

j=1;
for i = 1:size(shape1.WKS,2)
    shape1.fun_segment = [shape1.fun_segment repmat(ones(19248,1),1,size(shape1.parts,2)) .* shape1.parts];
    for k = 1:size(shape1.parts,2)
        W(j) = weightSegment;
        j = j+1; 
        W(j) = weightSegmentComm;
        j= j+1;
    end
    
end
for i = 1:size(shape1.HKS,2)
    shape1.fun_segment = [shape1.fun_segment repmat(ones(19248,1),1,size(shape1.parts,2)) .* shape1.parts];
    for k = 1:size(shape1.parts,2)
        W(j) = weightSegment;
        j = j+1; 
        W(j) = weightSegmentComm;
        j= j+1;
    end
end
for i = 1:size(shape2.WKS,2)
    shape2.fun_segment = [shape2.fun_segment repmat(shape2.WKS(:,i),1,size(shape2.parts,2)) .* shape2.parts];
end
for i = 1:size(shape2.HKS,2)
    shape2.fun_segment = [shape2.fun_segment repmat(shape2.HKS(:,i),1,size(shape2.parts,2)) .* shape2.parts];
end

for i = 1:size(shape1.HKS,2)
    W(j) = weightHKS;
    j = j+1;
    W(j) = weightHKSComm;
    j=j+1;
end

for i = 1:size(shape1.WKS,2)
    W(j) = weightWKS;
    j = j+1;
    W(j) = weightWKSComm;
    j=j+1;
end


shape1.fun = [shape1.fun_segment, shape1.HKS, shape1.WKS];
shape2.fun = [shape2.fun_segment, shape2.HKS, shape2.WKS];
%%

disp('Computing C');
C = calcCFromFuncsWeighted(shape1.fun,shape2.fun,diag(W),shape1.phi,shape2.phi,shape1.L,shape2.L);
disp('Done');
disp('Refining C');
%refinedC = refinementC(C,shape1,shape2,30);
refinedC = C;
disp('Refinement done');
disp('Computing point to point');


searchIndexParams = struct();
shape2toShape1 = flann_search(shape1.phi', refinedC'*shape2.phi', 1, searchIndexParams);
shape1toShape2 = flann_search(shape2.phi',C*shape1.phi',1,searchIndexParams);
disp('Done computing point to point');

%%
subplot(1,2,1);
options.face_vertex_color = shape1.vertex(shape2toShape1,2);
plot_mesh(shape2.vertex,shape2.faces,options);
subplot(1,2,2);
options2.face_vertex_color = shape1.vertex(:,2);
plot_mesh(shape1.vertex,shape1.faces,options2);

shading interp;
colormap jet;

%%

pointToHighlight = 15000;
fun = zeros(19248,1);
sigma = 10;
[D,S,Q] = perform_fast_marching_mesh(shape1.vertex,shape1.faces,pointToHighlight);

%%

subplot(1,2,1);
options.face_vertex_color = D(shape2toShape1);
plot_mesh(shape2.vertex,shape2.faces,options);
subplot(1,2,2);
options2.face_vertex_color = D;
plot_mesh(shape1.vertex,shape1.faces,options2);


shading interp;
colormap jet;

%%


geoError = zeros(19248,1);

for i = 1:19248
    options.end_points = shape2toShape1(i);
    [D,S,Q] = perform_fast_marching_mesh(shape1.vertex,shape1.faces,gt(i),options);
    err = D(shape2toShape1(i));
    geoError(i) = err;
    fprintf('%d \r', i/19248);
end

geoError2 = geoError;
for j = 0:0.01:40
   ind = int32(1 + j*100)
   pts(ind) = sum(geoError2<j);
end
plot(linspace(0,40,4001),pts/size(geoError,1))

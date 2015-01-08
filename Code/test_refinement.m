%% 2 NICO
%WKS + HKS + SEGMENT
%WKS + HKS + SEGMENT + REFINEMENT
%WKS + HKS + SEGMENT + REFINEMENT (N fois)
%% FOR ALL
% Graph de geodesic error *.mat file
% Images
% (+zoom intéressants)
% CHEVAL PAS GENTIL 1
% CHEVAL GENTIL 5
clear;
clc;
init;
mesh_name = '0003.isometry.5';
mesh_name2 = '0003.isometry.1';
name1 = 'Data/shrec10/0003.null.0.off';
name2 = strcat('Data/shrec10/',mesh_name,'.off');
gt_name = strcat('Data\shrec10gt\',mesh_name,'.labels');
gt = load(gt_name);
%%Avec la 4 c'est top horse
%%Avec la 5 c'est top horse
%%Avec la 2 c'est top dog
withRefinement = ones((100+1)*2,1);
withRefinement(101) = 0;
withRefinement(202) = 0;

%%
shape1 = getShape(name1);
shape2 = getShape(name2);
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
%%
shape1.indicComp(:,1) = 1.*(C1==11369);
shape2.indicComp(:,1) = 1.*(C2==15086);

nbIndicomp = 1;
for i = 1:size(list_matching,1)
    shape1.indicCompNOCONSTRAINTS(:,nbIndicomp) = 1.*(C1==list_matching(i,1));
    shape2.indicCompNOCONSTRAINTS(:,nbIndicomp) = 1.*(C2==list_matching(i,2));
    nbIndicomp = nbIndicomp +1
end

%%
% Test corresponding segment
disp('DO NOT FORGET TO CHANGE GT PATH');

shape1.parts = [];
shape2.parts = [];

shape1.parts = [shape1.parts 1.*(C1==15346)];
shape2.parts = [shape2.parts getAssociatedSegmentFromNull(1.*(C1==15346),gt)];
shape1.parts = [shape1.parts 1.*(C1==11369)];
shape2.parts = [shape2.parts getAssociatedSegmentFromNull(1.*(C1==11369),gt)];
shape1.parts = [shape1.parts 1.*(C1==4833)];
shape2.parts = [shape2.parts getAssociatedSegmentFromNull(1.*(C1==4833),gt)];

% figure();
% subplot(1,2,1);
% options1.face_vertex_color = segmentNull;
% plot_mesh(shape1.vertex,shape1.faces,options1);
% subplot(1,2,2);
% options2.face_vertex_color = segment4;
% plot_mesh(shape2.vertex,shape2.faces,options2);

clear shape1.fun_segment;
clear shape2.fun_segment;

shape1.fun_segment = [];
shape2.fun_segment = [];

W = zeros(1608,1);
weightWKS = 1;
weightHKS = 1;
weightSegment = 1;
weightWKSComm = 1;
weightHKSComm = 1;
weightSegmentComm = 1;

j=1;
for i = 1:size(shape1.WKS,2)
    shape1.fun_segment = [shape1.fun_segment repmat(shape1.WKS(:,i),1,size(shape1.parts,2)) .* shape1.parts];
    for k = 1:size(shape1.parts,2)
        W(j) = weightSegment;
        j = j+1;
        W(j) = weightSegmentComm;
        j= j+1;
    end
    
end
for i = 1:size(shape1.HKS,2)
    shape1.fun_segment = [shape1.fun_segment repmat(shape1.HKS(:,i),1,size(shape1.parts,2)) .* shape1.parts];
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
for iter = 1:numel(withRefinement)
    fprintf('performing %d/%d\n',iter,numel(withRefinement));
    if(iter == 102)
        mesh_name = mesh_name2;
        name2 = strcat('Data/shrec10/',mesh_name,'.off');
        gt_name = strcat('Data\shrec10gt\',mesh_name,'.labels');
        gt = load(gt_name);
        
        
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
        %%
        shape1.indicComp(:,1) = 1.*(C1==11369);
        shape2.indicComp(:,1) = 1.*(C2==15086);
        
        nbIndicomp = 1;
        for i = 1:size(list_matching,1)
            shape1.indicCompNOCONSTRAINTS(:,nbIndicomp) = 1.*(C1==list_matching(i,1));
            shape2.indicCompNOCONSTRAINTS(:,nbIndicomp) = 1.*(C2==list_matching(i,2));
            nbIndicomp = nbIndicomp +1
        end
        
        %%
        % Test corresponding segment
        disp('DO NOT FORGET TO CHANGE GT PATH');
        
        shape1.parts = [];
        shape2.parts = [];
        
        shape1.parts = [shape1.parts 1.*(C1==15346)];
        shape2.parts = [shape2.parts getAssociatedSegmentFromNull(1.*(C1==15346),gt)];
        shape1.parts = [shape1.parts 1.*(C1==11369)];
        shape2.parts = [shape2.parts getAssociatedSegmentFromNull(1.*(C1==11369),gt)];
        shape1.parts = [shape1.parts 1.*(C1==4833)];
        shape2.parts = [shape2.parts getAssociatedSegmentFromNull(1.*(C1==4833),gt)];
        
        % figure();
        % subplot(1,2,1);
        % options1.face_vertex_color = segmentNull;
        % plot_mesh(shape1.vertex,shape1.faces,options1);
        % subplot(1,2,2);
        % options2.face_vertex_color = segment4;
        % plot_mesh(shape2.vertex,shape2.faces,options2);
        
        clear shape1.fun_segment;
        clear shape2.fun_segment;
        
        shape1.fun_segment = [];
        shape2.fun_segment = [];
        
        W = zeros(1608,1);
        weightWKS = 1;
        weightHKS = 1;
        weightSegment = 1;
        weightWKSComm = 1;
        weightHKSComm = 1;
        weightSegmentComm = 1;
        
        j=1;
        for i = 1:size(shape1.WKS,2)
            shape1.fun_segment = [shape1.fun_segment repmat(shape1.WKS(:,i),1,size(shape1.parts,2)) .* shape1.parts];
            for k = 1:size(shape1.parts,2)
                W(j) = weightSegment;
                j = j+1;
                W(j) = weightSegmentComm;
                j= j+1;
            end
            
        end
        for i = 1:size(shape1.HKS,2)
            shape1.fun_segment = [shape1.fun_segment repmat(shape1.HKS(:,i),1,size(shape1.parts,2)) .* shape1.parts];
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
    end
    
    %%
    disp('Computing C');
    C = calcCFromFuncsWeighted(shape1.fun,shape2.fun,diag(W),shape1.phi,shape2.phi,shape1.L,shape2.L);
    disp('Done');
    if( withRefinement(iter) )
        disp('Refining C');
        refinedC = refinementC(C,shape1,shape2,30);
        disp('Refinement done');
    end
    
    disp('Computing point to point');
    searchIndexParams = struct();
    shape2toShape1 = flann_search(shape1.phi', refinedC'*shape2.phi', 1, searchIndexParams);
    shape1toShape2 = flann_search(shape2.phi',C*shape1.phi',1,searchIndexParams);
    disp('Done computing point to point');
    
    
    
    %%
    % pointToHighlight = 15000;
    % fun = zeros(19248,1);
    % sigma = 10;
    % [D,S,Q] = perform_fast_marching_mesh(shape1.vertex,shape1.faces,pointToHighlight);
    
    
    
    geoError = zeros(19248,1);
    disp('Performing geodesic error');
    tic
    for i = 1:19248
        options.end_points = shape2toShape1(i);
        [D,S,Q] = perform_fast_marching_mesh(shape1.vertex,shape1.faces,gt(i),options);
        err = D(shape2toShape1(i));
        geoError(i) = err;
        fprintf('%d \r', i/19248);
    end
    toc
    %%
    filename = strcat('./test/2_Refinement/geo_error.',mesh_name,'.mat');
    save(filename,'geoError');
end
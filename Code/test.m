%clear; clc;
init;

name1 = 'Data/shrec10/0002.null.0.off';
name2 = 'Data/shrec10/0002.isometry.1.off';

shape1 = getShape(name1);
shape2 = getShape(name2);

%%
%Compute connected component of shape 1

C1 = persistance_based_segmentation(shape1,7);
shape1.connected_component = C1;
list_label_C1 = union(C1,C1);

%%
%Compute connected component of shape 2

C2 = persistance_based_segmentation(shape2,7);
shape2.connected_component = C2;
list_label_C2 = union(C2,C2);

%%
%Matching segment bewteen shape1 and shape2
list_matching = getSegmentMatching(shape1,shape2);

%%
%Visually check that each segment correspondence is good
for i = 1:size(list_matching,1)
    fprintf('Segments matching num %d with err=%f\n',i,list_matching(i,3));
    figure(1);
    options.face_vertex_color = compute_color_from_connected_component(C1, list_matching(list_matching(:,1)==list_matching(i,1),1));
    plot_mesh(shape1.vertex,shape1.faces,options);
    shading interp; colormap jet(256);
    
    figure(2);
    options.face_vertex_color = compute_color_from_connected_component(C2, list_matching(list_matching(:,2)==list_matching(i,2),2));
    plot_mesh(shape2.vertex,shape2.faces,options);
    shading interp; colormap jet(256);
    fprintf('Press start too check next correspondence\n');
    pause();
end

%% 

nbIndicomp = 1;
for i = 1:size(list_matching,1)
    fprintf('Couple %d, C1 : %d C2 : %d\n',i,sum(C1==list_matching(i,1)),sum(C2==list_matching(i,2)));
    diff_num_vertex = abs(sum(C1==list_matching(i,1))-sum(C2==list_matching(i,2)))/(sum(C1==list_matching(i,1))+sum(C2==list_matching(i,2)));
    fprintf('diff of nb vertex normalized : %f\n',diff_num_vertex);
    if( list_matching(i,4) < 0.1 && sum(C1==list_matching(i,1)) > 150)
        shape1.indicComp(:,i) = 1.*(C1==list_matching(i,1));
        shape2.indicComp(:,i) = 1.*(C2==list_matching(i,2));
        nbIndicomp = nbIndicomp +1
    end
end


%%
[C] = getFunctionalMap(shape1,shape2);


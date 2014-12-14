function list_matching = getSegmentMatching( shape1, shape2 )

list_matching = [];

C1 = shape1.connected_component;
C2 = shape2.connected_component;
list_label_C1 = union(C1,C1);
list_label_C2 = union(C2,C2);

list_descriptors_C1 = compute_descriptors_for_matching(shape1);
list_descriptors_C2 = compute_descriptors_for_matching(shape2);

%[~,perm1] = sort(list_descriptors_C1);
%[~,perm2] = sort(list_descriptors_C2);
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
    diff_nb_vertex = abs(sum(C1==list_label_C1(idx_row(idx_col)))-sum(C2==list_label_C2(idx_col)))...
        /(sum(C1==list_label_C1(idx_row(idx_col)))+sum(C2==list_label_C2(idx_col)));
    list_matching(i,:) = [ list_label_C1(idx_row(idx_col)) list_label_C2(idx_col) diff diff_nb_vertex ];
    
    % Replace idx_row(idx_col)th row and idx_colth col of err with a high
    % value
    err(idx_row(idx_col),:) = INFTY;
    err(:,idx_col) = INFTY;
end

end


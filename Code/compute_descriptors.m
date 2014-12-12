function list_descriptors_C = compute_descriptors(shape)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
C = shape.connected_component;

list_label_C = union(C,C);
list_descriptors_C = zeros(size(list_label_C,1),1);

% WKS with a fixed energy e=5
log_E=log(max(abs(shape.E),1e-6))';
e=linspace(log_E(2),(max(log_E))/1.02,size(shape.WKS,2));
[~,idx] = min(abs(e-5));

for i=1:size(list_label_C,1)
    list_descriptors_C(i) = sum(shape.WKS(C==list_label_C(i),idx));
end
end


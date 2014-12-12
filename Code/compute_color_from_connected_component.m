function color = compute_color_from_connected_component( C )
color = zeros(size(C,1),1);
component_value = union(C,C);
for i=1:size(component_value,1)
    color(C==component_value(i)) = i*2500;
end
end


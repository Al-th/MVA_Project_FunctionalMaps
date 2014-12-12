function color = compute_color_from_connected_component( C , component_value )
color = zeros(size(C,1),1);
for i=1:size(component_value,1)
    color(C==component_value(i)) = i*2000;
end
end


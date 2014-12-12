function connected_component = persistance_based_segmentation( shape, nb_wanted_seg )
%This function segment the shape in maximum nb_wanted_seg parts using the
%algorithm describe in [Skraba et al 2010] 
%
% Find the correct tau thanks to the persistence diagram PD
fprintf('Computing Persistence Diagram to determined tau for %d components...',nb_wanted_seg);
[~,PD] = compute_connected_component(shape,100000);
diff = sort((PD(:,1)-PD(:,2)),'descend');
nb_wanted_seg = nb_wanted_seg -1;
if(numel(diff) > nb_wanted_seg)
    while(nb_wanted_seg > 1 & diff(nb_wanted_seg) < mean(diff))
        nb_wanted_seg = nb_wanted_seg - 1;
    end
end

tau = diff(nb_wanted_seg)-1e-6;
fprintf('done. \n');

fprintf('Computing connected component with tau = %f...',tau);
[connected_component,~] = compute_connected_component(shape,tau);
fprintf('done. \n');
end


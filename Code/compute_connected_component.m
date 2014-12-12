function [connected_component, PD] = compute_connected_component(shape,tau)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%parameters

num_vertices = size(shape.vertex,1);
num_faces = size(shape.faces,1);
C = zeros(num_vertices,1);
PD = [];
nb_merge = 0;
% sort vertex w.r.t to WKS
wks = shape.WKS(:,5); %sum(shape.WKS');
[~,perm]=sort(wks,'descend');
perm = perm'; % shape.vertex(perm(1))

for i=1:num_vertices,
    %process vertex x
    idx_x = perm(i);
    x = shape.vertex(idx_x,:);
    
    %%
    %Compute list of vertex in a one-ring neighborhood of x
    
    %find faces containing x
    idx_faces = [];
    for k=1:3
        cur_idx = find(shape.faces(:,k)==idx_x);
        idx_faces = [idx_faces cur_idx'];
    end
    % Recupere les vertex de chaque face et supprime les doublons
    tmp = shape.faces(idx_faces,:);
    ngb_x_index = intersect(tmp(tmp~=idx_x),tmp(tmp~=idx_x));
    
    %%
    %Determine if x is a local maximum w.r.t wks
    local_max = wks(idx_x);
    idx_local_max = idx_x;
    for k=ngb_x_index',
        if( wks(k) > local_max)
            local_max = wks(k);
            idx_local_max = k;
        end
    end
    if(idx_local_max == idx_x)
        % Create a new component
        C(idx_x) = idx_x;
    else
        %Add x to the component with the highest wks value
        C(idx_x) = C(idx_local_max);
    end
    
    %%
    %If the vertex x is adjacent to two or more existing components,
    %we check the persistence of the components and merge them
    %only if they are not tau-persistent
    nb_adj_component = numel(find(C(ngb_x_index)));
    a = wks(idx_x);
    if( nb_adj_component > 1 );
        
        x2 = wks(C(idx_x));
        for k = ngb_x_index'
            if(C(k) ~= 0)
                x1 = wks(C(k));
                if( x1 < x2 )
                    if( x1 - a <= tau )
                        %We merge
                        nb_merge = nb_merge +1;
                        C(C==C(k)) = C(idx_x);
                        PD(nb_merge,:) = [x1 a];
                    else
                        %We don't merge
                    end
                end
            end
        end
    end
end

connected_component = C;

end


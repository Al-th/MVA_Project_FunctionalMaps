clear; clc;
init;

%%
%load shape 1
name = 'Data/shrec10/0003.null.0.off';
shape1.name = name;
[vertex,faces] = read_off(name);
vertex = vertex'; faces = faces';
shape1.vertex = vertex;
shape1.faces = faces;
%Compute ALB
[PHI,E,L,Am] = ALB_spectrum(vertex,faces);
shape1.phi = PHI;
shape1.E = E;
shape1.Am = Am;
%Compute HKS
hks = HKS(PHI, E, diag(Am), false);
shape1.HKS = hks;
%Compute WKS
WKS = compute_WKS_from_EV(E,PHI);
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
[PHI,E,L,Am] = ALB_spectrum(vertex,faces);
shape2.phi = PHI;
shape2.E = E;
shape2.Am = Am;
%Compute HKS
hks = HKS(PHI, E, diag(Am), false);
shape2.HKS = hks;
%Compute WKS
WKS = compute_WKS_from_EV(E,PHI);
shape2.WKS = WKS;
clear WKS; clear hks;
clear PHI; clear E; clear L; clear Am;
clear vertex; clear faces; clear name;


%%
%Compute persistence-based segmentation based on [Skraba et al. 2010]

%parameters
tau = 5;
num_vertices = size(vertices,1);
num_faces = size(faces,1);
C = zeros(num_vertices);

% sort vertex w.r.t to WKS
wks = sum(shape1.WKS');
[~,perm]=sort(wks,'descend');
perm = perm'; % shape1.vertex(perm(1))

for i=1:num_vertices,
   %process vertex x 
    idx_x = perm(i);
    x = shape1.vertex(idx_x,:);
    
    %%
    %Compute list of vertex in a one-ring neighborhood of x
    
    %find faces containing x
    idx_faces = [];
    for k=1:3      
       cur_idx = find(shape1.faces(:,k)==idx_x);
       idx_faces = [idx_faces cur_idx'];       
    end
    % Recupere les vertex de chaque face et supprime les doublons
    tmp = shape1.faces(idx_faces,:);
    ngb_x_index = intersect(tmp(tmp~=idx_x),tmp(tmp~=idx_x));
    
    %%
    %Determine if x is a local maximum w.r.t wks
    local_max = wks(idx_x);
    idx_local_max = idx_x;
    for k=ngb_x_index,
        if( wks(k) > local_max)
            local_max = wks(k);
            idx_local_max = k;
        end
    end
    if(idx_local_max == idx_x)
        % Create a new component
        C(id_x) = idx_x;
    else
        %Add x to the component with the highest wks value
        C(id_x) = C(idx_local_max);
    end
    
    %%
    %If the vertex x is adjacent to two or more existing components, 
    %we check the persistence of the components and merge them 
    %only if they are not tau-persistent
    nb_adj_component = numel(find(C(ngb_x_index)));
    a = wks(idx_x);
    if( nb_adj_component > 1 );
        
        x2 = wks(C(idx_x));
        
        for k = ngb_x_index
            x1 = wks(C(k));
            if( x1 < x2 )
                if( x1 - a <= tau )
                    %We merge
                    C(C==C(k)) = C(idx_x);
                else
                    %We don't merge
                end
            end
        end
    end
end





%%
% Add constraint to linear system
a1 = shape1.phi'*shape1.Am*shape1.HKS;
a2 = shape1.wks_phi'*shape1.WKS;
a = [a1 a2];

b1 = shape2.phi'*shape2.Am*shape2.HKS;
b2 = shape2.wks_phi'*shape2.WKS;
b = [b1 b2];




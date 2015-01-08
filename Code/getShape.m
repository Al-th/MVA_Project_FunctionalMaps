function shape = getShape(name,alpha)

alboptions.n_eigenvalues = 40; 
if nargin > 1
   alboptions.alpha = alpha 
end
shape.nbEigenvalues = alboptions.n_eigenvalues;    
shape.name = name;
[vertex,faces] = read_off(shape.name);
vertex = vertex'; faces = faces';
shape.vertex = vertex;
shape.faces = faces;
%Compute ALB
[PHI,E,L,Am] = ALB_spectrum(vertex,faces,alboptions);
shape.phi = PHI;
shape.E = E;
shape.Am = Am;
shape.L  = L;
%Compute HKS
hks = HKS(PHI, E, diag(Am), false);
shape.HKS = hks;
%Compute WKS
WKS = compute_WKS_from_EV(E,PHI,alboptions);
shape.WKS = WKS;


end
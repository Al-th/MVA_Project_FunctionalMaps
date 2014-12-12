function [constraints,equalTo] = computeOperatorCommutativityConstraints(f,shape1,shape2)
%CA = B*C*D

A = shape1.phi'*shape1.Am*shape1.L*f;
B = shape2.phi'*shape2.Am*shape2.L*shape2.phi;
D = shape1.phi'*shape1.Am*f;

n = size(shape1.phi,2);
MA = zeros(n,n*n);
for i = 1:n
    for j = 1:n
        MA(i,(i-1)*n+j) = A(j); 
    end
end

MBD = zeros(n,n*n);

for i = 1:n
    for j = 1:n*n
        dIndex = mod(j-1,n)+1;
        M1(i,j) = -D(dIndex) * B(i,floor((j-1)/n)+1); 
    end
end
constraints = MA+MBD;

equalTo = zeros(n,1);

end
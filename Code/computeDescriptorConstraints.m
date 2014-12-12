function [constraints, equalTo] = computeDescriptorConstraints(f1,shape1,f2,shape2)

A = shape1.phi'*shape1.Am*f1;
B = shape2.phi'*shape2.Am*f2;

n = size(shape1.phi,2);
MA = zeros(n,n*n);
for i = 1:n
    for j = 1:n
        MA(i,(i-1)*n+j) = A(j); 
    end
end
constraints = MA;
equalTo = B;
end
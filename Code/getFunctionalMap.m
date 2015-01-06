function [C] = getFunctionalMap(shape1,shape2)
    %%Writing linear constraints for C determination

    A = [];
    b = [];
    
%     %Adding phi constraints
%     for i = 1:size(shape1.phi,2)
%        [dA,db] = computeDescriptorConstraints(shape1.phi(:,i),shape1, shape2.phi(:,i),shape2);
%         A = sparse([A ; dA]);
%         b = sparse([b ; db]);
%     end
%     

    
    tic
    %%add descriptor preservation constraints
    for i = 1:size(shape1.fun,2)
        [dA,db] = computeDescriptorConstraints(shape1.fun(:,i),shape1,shape2.fun(:,i),shape2);
        A = sparse([A ; dA]);
        b = sparse([b ; db]);
        
%         [dA,db] = computeOperatorCommutativityConstraints(shape1.fun(:,i),shape1,shape2);
%         A = sparse([A; dA]);
%         b = sparse([b; db]);
    end
    toc


    disp('Solving linear system Ax=b')
    tic
    x = lsqlin(A,b);
    %x = lsqr(A,b);
    %x = mldivide(A,b);
    toc
    
    n = shape1.nbEigenvalues;
    C = reshape(x,n,n)';
end
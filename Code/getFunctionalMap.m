function [C] = getFunctionalMap(shape1,shape2)
    %%Writing linear constraints for C determination

    A = [];
    b = [];


    disp('Adding HKS preservation constraints')
    tic
    %%add descriptor preservation constraints
    for i = 1:size(shape1.HKS,2)
        [dA,db] = computeDescriptorConstraints(shape1.HKS(:,i),shape1,shape2.HKS(:,i),shape2);
        A = sparse([A ; dA]);
        b = sparse([b ; db]);
    end
    toc

    disp('Adding WKS preservation constraints')
    tic
    for i = 1:size(shape1.WKS,2)
        [dA,db] = computeDescriptorConstraints(shape1.WKS(:,i),shape1,shape2.WKS(:,i),shape2);
        A = sparse([A ; dA]);
        b = sparse([b ; db]);
    end
    toc
    
    disp('Adding segment preservation constraints')
    tic
    for i = 1:size(shape1.indicComp,2)
        [dA,db] = computeDescriptorConstraints(shape1.indicComp(:,i),shape1,shape2.indicComp(:,i),shape2);
        A = sparse([A ; dA]);
        b = sparse([b ; db]);
    end
    toc


    disp('Adding HKS operator commutativity constraints')
    tic
    %%add operator commutativity constraints
    for i = 1:size(shape1.HKS,2)
        [dA,db] = computeOperatorCommutativityConstraints(shape1.HKS(:,i),shape1,shape2);
        A = sparse([A; dA]);
        b = sparse([b; db]);
    end
    toc

    disp('Adding WKS operator commutativity constraints')
    tic
    for i = 1:size(shape1.WKS,2)
        [dA,db] = computeOperatorCommutativityConstraints(shape1.WKS(:,i),shape1,shape2);
        A = sparse([A; dA]);
        b = sparse([b; db]);
    end
    toc
    
    disp('Adding segment operator commutativity constraints')
    tic
    for i = 1:size(shape1.indicComp,2)
        [dA,db] = computeOperatorCommutativityConstraints(shape1.indicComp(:,i),shape1,shape2);
        A = sparse([A; dA]);
        b = sparse([b; db]);
    end
    toc


    %
    disp('Solving linear system Ax=b')
    tic
    %x = lsqlin(A,b);
    %x = lsqr(A,b);
    x = mldivide(A,b);
    toc
    
    n = shape1.nbEigenvalues;
    C = reshape(x,n,n)';
end
function [C] = getFunctionalMapWeighted(shape1,shape2)
    %%Writing linear constraints for C determination

    A = [];
    b = [];
    w = [];
    
    %%Adding phi constraints
    weight = 1;
    for i = 1:size(shape1.phi,2)
       [dA,db] = computeDescriptorConstraints(shape1.phi(:,i),shape1, shape2.phi(:,i),shape2);
        A = sparse([A ; dA]);
        b = sparse([b ; db]);
        w = [w; weight.*ones(size(db,1),1)];
    end

    disp('Adding HKS preservation constraints')
    weight = 0.1;
    tic
    %%add descriptor preservation constraints
    for i = 1:size(shape1.HKS,2)
        [dA,db] = computeDescriptorConstraints(shape1.HKS(:,i),shape1,shape2.HKS(:,i),shape2);
        A = sparse([A ; dA]);
        b = sparse([b ; db]);
        w = [w; weight.*ones(size(db,1),1)];
    end
    toc

    disp('Adding WKS preservation constraints')
    weight = 0.1;
    tic
    for i = 1:size(shape1.WKS,2)
        [dA,db] = computeDescriptorConstraints(shape1.WKS(:,i),shape1,shape2.WKS(:,i),shape2);
        A = sparse([A ; dA]);
        b = sparse([b ; db]);
        w = [w; weight.*ones(size(db,1),1)];
    end
    toc
    
    disp('Adding segment preservation constraints')
    weight = 1.;
    tic
    for i = 1:size(shape1.indicComp,2)
        [dA,db] = computeDescriptorConstraints(shape1.indicComp(:,i),shape1,shape2.indicComp(:,i),shape2);
        A = sparse([A ; dA]);
        b = sparse([b ; db]);
        w = [w; weight.*ones(size(db,1),1)];
    end
    toc


    disp('Adding HKS operator commutativity constraints')
    weight = 0.1;
    tic
    %%add operator commutativity constraints
    for i = 1:size(shape1.HKS,2)
        [dA,db] = computeOperatorCommutativityConstraints(shape1.HKS(:,i),shape1,shape2);
        A = sparse([A; dA]);
        b = sparse([b; db]);
        w = [w; weight.*ones(size(db,1),1)];
    end
    toc

    disp('Adding WKS operator commutativity constraints')
    weight = 0.1;
    tic
    for i = 1:size(shape1.WKS,2)
        [dA,db] = computeOperatorCommutativityConstraints(shape1.WKS(:,i),shape1,shape2);
        A = sparse([A; dA]);
        b = sparse([b; db]);
        w = [w; weight.*ones(size(db,1),1)];
    end
    toc
    
    disp('Adding segment operator commutativity constraints')
    weight = 1.;
    tic
    for i = 1:size(shape1.indicComp,2)
        [dA,db] = computeOperatorCommutativityConstraints(shape1.indicComp(:,i),shape1,shape2);
        A = sparse([A; dA]);
        b = sparse([b; db]);
        w = [w; weight.*ones(size(db,1),1)];
    end
    toc
    
    W = diag(w);

    %
    disp('Solving linear system Ax=b')
    tic
    %x = lsqlin(A,b);
    %x = lsqr(A,b);
    x = mldivide(A,b);
    toc
    
    %
    disp('Solving linear system AtWtWAx=AtWtWb')
    tic
    %x = lsqlin(A,b);
    %x = lsqr(A,b);
    x1 = mldivide(A'*(W'*(W*A)),A'*(W'*(W*b)));
    toc
    
    
    n = shape1.nbEigenvalues;
    C = reshape(x1,n,n)';
end
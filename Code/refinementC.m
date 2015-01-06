function [ refinedC ] = refinementC( C, shape1, shape2, numRefinements )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

C = C';
refinedC = C;
prevC = C;
b1 = shape1.phi';
b2 = shape2.phi';
prevVal = Inf;
searchIndexParams = struct();
shape1toshape2 = [];
timeit.ICP = 0;
for icpIter = 0:numRefinements
    
    b2Perm = refinedC * b1;
    b1Perm = refinedC'* b2;
    if icpIter == numRefinements
        %searchIndexParams = struct('algorithm', 'linear');
    else
        searchIndexParams = struct();
    end
    
    tic;
    %
    prevShape1toshape2 = shape1toshape2;
    shape1toshape2 = flann_search(b2, b2Perm, 1, searchIndexParams);
    %shape2toshape1 = flann_search(refinedC * basis1', b2, 1, searchIndexParams);
    B = b2(:, shape1toshape2);
    A = b1;
    newVal = norm(refinedC * A - B,'fro');
    if prevVal < newVal
        fprintf('   Refine worsened the results restoring prevC\n');
        refinedC = prevC;
        shape1toshape2 = prevShape1toshape2;
        b2Perm = refinedC * b1;
        b1Perm = refinedC'* b2;
        break;
    end
    prevVal = newVal;
    
    flannTime = toc;
    timeit.ICP = timeit.ICP + flannTime;
    %fprintf('  flann time = %f secs', flannTime);
    if icpIter < numRefinements
        fprintf('  Preforming refinement (%d/%d)\n', icpIter + 1, numRefinements);
        
        tic;
        %find refinedC = argmin ||refinedC * basis1 - basis2(:, shape1toshape2result)||
        % s.t. refinedC * refinedC' = I
        
        [U,S,V] = svd(B*A');
        U(:,end)=U(:,end)*det(U*V');
        newCR=U*V';
        timeit.ICP = timeit.ICP + toc;
        
        
        crdiff = sum(sum((abs(newCR - refinedC))));
        fprintf('    crdiff = %f\n', crdiff);
        
        if crdiff < 0.01
            fprintf('  Refinments converged stopping refinments\n')
            break;
        end
        prevC = refinedC;
        refinedC = newCR;
    end
end

end


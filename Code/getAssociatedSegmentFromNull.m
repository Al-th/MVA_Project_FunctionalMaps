function segment = getAssociatedSegment(originSegment, groundTruth)

segment = zeros(size(originSegment));
for i = 1:size(originSegment,1)
    if originSegment(i) == 1
       %Si la vertex appartient à un segment, on cherche le correspondant
       %dans la GT
       segment(find(groundTruth == i)) = 1;
    end
end

end
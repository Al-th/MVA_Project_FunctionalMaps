function distance = computeDistanceFrom(vertex,a,b,c)
    distance = sqrt((vertex(:,1)-a).^2+(vertex(:,2)-b).^2+(vertex(:,3)-c).^2);
end
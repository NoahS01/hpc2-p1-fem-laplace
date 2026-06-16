function kFound = findTriangleIndexLocal( ...
    p, coordinates, elements, node2elem, rSearch, localNodes)

%% Nahe Knoten suchen

dist = vecnorm(coordinates(localNodes,:) - p, 2, 2);
nearLocal = localNodes(dist <= rSearch);

if isempty(nearLocal)
    error('No nodes found in local search radius');
end

%% Kandidatendreiecke sammeln

candidates = unique([node2elem{nearLocal}]);

%% Kandidaten testen

for i = 1:length(candidates)

    k = candidates(i);

    Pk = coordinates(elements(k,:),:);

    if pointInTriangle(p, Pk)
        kFound = k;
        return;
    end

end

error('Point not found in local triangle search');

end
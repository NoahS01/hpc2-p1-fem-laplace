function [integral_u, integral_dn] = computeCircleIntegrals( ...
    coordinates, elements, u, cx, cy, radius, N, hMesh)

%% Kreis diskretisieren

theta = linspace(0, 2*pi, N+1);
theta = theta(1:N);
circle = [cx + radius*cos(theta); cy + radius*sin(theta)]';

ds = 2*pi*radius / (N);

%% Suchparameter

rSearch = 2*hMesh;

%% Suchstruktur vorbereiten

node2elem = buildNode2Elem(elements, size(coordinates,1));

xmin = cx - radius - rSearch;
xmax = cx + radius + rSearch;
ymin = cy - radius - rSearch;
ymax = cy + radius + rSearch;

localNodes = find( ...
    coordinates(:,1) >= xmin & coordinates(:,1) <= xmax & ...
    coordinates(:,2) >= ymin & coordinates(:,2) <= ymax);

%% Elementkoeffizienten vorberechnen

coeffElem = computeElementCoefficients(coordinates, elements, u);

%% Kreiswerte berechnen

kCircle = zeros(N,1);
kLast = 0;

for i = 1:N
    p = circle(i,:);

    k = 0;

    % zuerst letztes Dreieck testen
    if kLast > 0
        PLast = coordinates(elements(kLast,:),:);

        if pointInTriangle(p, PLast)
            k = kLast;
        end
    end

    % lokale Suche als Fallback
    if k == 0
        k = findTriangleIndexLocal( ...
            p, coordinates, elements, node2elem, rSearch, localNodes);
    end

    kCircle(i) = k;
    kLast = k;
end

% Update: vektorisierte Berechnung, statt in for-Schleife
coeff = coeffElem(kCircle,:);

% Funktionswerte
uCircle = coeff(:,1) ...
        + coeff(:,2).*circle(:,1) ...
        + coeff(:,3).*circle(:,2);

% Normalenvektoren auf dem Kreis
normals = [(circle(:,1)-cx)/radius, ...
           (circle(:,2)-cy)/radius];

% Normalenableitung
dnCircle = coeff(:,2).*normals(:,1) ...
         + coeff(:,3).*normals(:,2);



%% Integral von u über den Kreis

integral_u  = ds * sum(uCircle);

%% Integral der Normalenableitung

integral_dn = ds * sum(dnCircle);

end

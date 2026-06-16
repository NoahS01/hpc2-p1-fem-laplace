function [integral_u, integral_dn] = computeCircleIntegrals( ...
    coordinates, elements, u, cx, cy, radius, N, hMesh)

%% Kreis diskretisieren

theta = linspace(0, 2*pi, N);
circle = [cx + radius*cos(theta); cy + radius*sin(theta)]';

ds = 2*pi*radius / (N-1);

%% Suchparameter

rSearch = 2*hMesh;
buffer = 2*hMesh;

%% Suchstruktur vorbereiten

node2elem = buildNode2Elem(elements, size(coordinates,1));

xmin = cx - radius - buffer;
xmax = cx + radius + buffer;
ymin = cy - radius - buffer;
ymax = cy + radius + buffer;

localNodes = find( ...
    coordinates(:,1) >= xmin & coordinates(:,1) <= xmax & ...
    coordinates(:,2) >= ymin & coordinates(:,2) <= ymax);

%% Elementkoeffizienten vorberechnen

coeffElem = computeElementCoefficients(coordinates, elements, u);

%% Kreiswerte berechnen

uCircle = zeros(N,1);
dnCircle = zeros(N,1);

kLast = 0;

for i = 1:N

    p = circle(i,:);

    k = 0;

    % Zuerst prüfen, ob der Punkt noch im zuletzt gefundenen Dreieck liegt
    if kLast > 0

        PLast = coordinates(elements(kLast,:),:);

        if pointInTriangle(p, PLast)
            k = kLast;
        end

    end

    % Falls nicht: lokale Suche
    if k == 0

        k = findTriangleIndexLocal( ...
            p, coordinates, elements, node2elem, rSearch, localNodes);

    end

    kLast = k;

    coeff = coeffElem(k,:);

    % u_h(x,y) = a + b*x + c*y
    uCircle(i) = coeff(1) + coeff(2)*p(1) + coeff(3)*p(2);

    % Normalenvektor radial nach außen
    n_vec = (p - [cx, cy]) / radius;

    % grad(u_h) = [b, c]
    dnCircle(i) = coeff(2)*n_vec(1) + coeff(3)*n_vec(2);

end

%% Integral von u über den Kreis

sum_u = 0;

for i = 1:N-1
    sum_u = sum_u + 0.5 * (uCircle(i) + uCircle(i+1));
end

integral_u = sum_u * ds;

%% Integral der Normalenableitung

sum_dn = sum(dnCircle(1:N-1));

integral_dn = sum_dn * ds;

end
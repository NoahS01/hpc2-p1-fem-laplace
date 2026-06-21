%% Geometrie und Mesh

hMesh = 0.3;

tStart = tic;
tMesh = tic;
[coordinates, elements, boundary, dirichlet, neumann] = createMesh(false, hMesh);
tMesh = toc(tMesh);

%% Kreis für Linienintegrale

radius = 2;
cx = 2.5;
cy = -1.5;
N = 2000;

%% Randbedingungen

f = @(x) zeros(size(x,1),1);
g = @(x) zeros(size(x,1),1);

nodes5 = unique(boundary{6}(:));
nodes3 = unique(boundary{13}(:));
nodes7 = unique([boundary{14}(:); boundary{15}(:)]);

coords5 = coordinates(nodes5,:);
coords3 = coordinates(nodes3,:);
coords7 = coordinates(nodes7,:);

uD = @(x) uDirichlet(x, coords5, coords3, coords7);

%% PDE lösen

tSolve = tic;
[x, Energie] = solveLaplace(coordinates, elements, dirichlet, neumann, f, g, uD);
tSolve = toc(tSolve);

%% Linienintegrale berechnen

tPost = tic;
[integral_u, integral_dn] = computeCircleIntegrals( ...
    coordinates, elements, x, cx, cy, radius, N, hMesh);
tPost = toc(tPost);

tTotal = toc(tStart);

%% Laufzeit-Ausgabe

fprintf('Anzahl Knoten:          %d\n', size(coordinates,1));
fprintf('Anzahl Elemente:        %d\n', size(elements,1));
fprintf('Netzerzeugung:          %.4f s\n', tMesh);
fprintf('FEM-Assembly + Lösung:  %.4f s\n', tSolve);
fprintf('Postprocessing:         %.4f s\n', tPost);
fprintf('Gesamtlaufzeit:         %.4f s\n', tTotal);
fprintf('Anteil Postprocessing:  %.1f %%\n', 100*tPost/tTotal);
clc
clear
close all

%% Geometrie und Mesh

hMesh = 0.3;
[coordinates, elements, boundary, dirichlet, neumann] = createMesh(true, hMesh);

%% Kreis für Linienintegrale

radius = 2;
cx = 2.5;
cy = -1.5;
N = 2000;

%% Randbedingungen

f = @(x) zeros(size(x,1),1);
g = @(x) zeros(size(x,1),1);

nodes5 = unique(boundary{13}(:));
nodes3 = unique(boundary{6}(:));
nodes7 = unique([boundary{14}(:); boundary{15}(:)]);

coords5 = coordinates(nodes5,:);
coords3 = coordinates(nodes3,:);
coords7 = coordinates(nodes7,:);

uD = @(x) uDirichlet(x, coords5, coords3, coords7);

%% PDE lösen

[x, Energie] = solveLaplace(coordinates, elements, dirichlet, neumann, f, g, uD);

%% Linienintegrale berechnen

[integral_u, integral_dn] = computeCircleIntegrals( ...
    coordinates, elements, x, cx, cy, radius, N, hMesh);

%% Plot der FEM-Lösung

theta = linspace(0, 2*pi, N);
circle = [cx + radius*cos(theta); cy + radius*sin(theta)]';

figure;
trisurf(elements, coordinates(:,1), coordinates(:,2), x);

axis equal;
view(3);
shading interp;
colorbar;

title('P1-FEM Lösung der Laplace-Gleichung');
xlabel('x');
ylabel('y');
zlabel('u_h');

hold on
plot(circle(:,1), circle(:,2), 'r', 'LineWidth', 2)
hold off

%% Ausgabe

disp(['Energie: ', num2str(Energie)]);
disp(['I_u: ', num2str(integral_u)]);
disp(['I_dn: ', num2str(integral_dn)]);

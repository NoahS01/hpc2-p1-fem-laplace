clc
clear
close all


%% Schritt- / Gitterweite
mesh_sizes          = [2:-0.05:0.20, 0.18:-0.02:0.02 ];
circle_stepsizes    = [6: -0.1: 1.1,    1:-0.02:0.02 ];

%% Ergebnisse / Daten
num_of_Nodes    = zeros(length(mesh_sizes), 1);
num_of_C_Nodes  = zeros(length(circle_stepsizes), 1);
Energies        = zeros(length(mesh_sizes), 1);

circle_u        = [];
integrals_u     = zeros(length(mesh_sizes), length(circle_stepsizes));
integrals_dudn  = zeros(length(mesh_sizes), length(circle_stepsizes));

% stop_watch
sw_meshing          = zeros(length(mesh_sizes), 1); % Dauer des meshings
sw_solvePDE         = zeros(length(mesh_sizes), 1); % Dauer des lösens
sw_solveIntegrals   = zeros(length(mesh_sizes), length(circle_stepsizes));

%% Iteration über die Gitterweite des Traktors
for idx_mesh = 1:length(mesh_sizes)
    hMesh = mesh_sizes(idx_mesh)

    %% Traktor Geometrie
    tic
    [coordinates,elements,boundary,dirichlet,neumann] = createMesh(false, hMesh);
    sw_meshing(idx_mesh) = toc;
    num_of_Nodes(idx_mesh) = size(coordinates, 1);           

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
    tic
    [x,Energie] = solveLaplace( ...
        coordinates, elements, dirichlet, neumann, f, g, uD);
    sw_solvePDE(idx_mesh) = toc;

    Energies(idx_mesh) = Energie;

    %% Iteration über den diskretisierten Kreis
    for idx_circle = 1:length(circle_stepsizes)
        h_circle = circle_stepsizes(idx_circle);

        %% Kreis Geometrie
        radius = 2.0;
        cx = 2.5; cy = -1.5;
        N = round(2*pi*radius/h_circle);
        num_of_C_Nodes(idx_circle) = N;

        %% Linienintegrale berechnen   
        tic
        [integral_u, integral_dudn] = computeCircleIntegrals( ...
            coordinates, elements, x, cx, cy, radius, N, hMesh );
        sw_solveIntegrals(idx_mesh, idx_circle) = toc;
        %circle_u(end+1,:)                       = uCircle;
        integrals_u(idx_mesh, idx_circle)       = integral_u;
        integrals_dudn(idx_mesh, idx_circle)    = integral_dudn;

    end
end


%% Speichern in Datei

writematrix(integrals_u,'csv/convergence_Iu.csv')
writematrix(integrals_dudn,'csv/convergence_Idudn.csv')

%% Plot der Ergebnisse

% Zielordner für Graphen
graphDir = './Graphs';
if ~exist(graphDir, 'dir')
    mkdir(graphDir);
end

[X,Y] = meshgrid(num_of_Nodes, num_of_C_Nodes);

% Konvergenzplot für FEM mit Energie über Energie(h_max)
fig = figure('Units','centimeters','Position',[5 5 24 12]);
plot(num_of_Nodes, Energies, 'LineWidth', 1.5)
xlabel("DoF Mesh")
ylabel("Energie")
set(gca,'XScale','log')
exportgraphics(fig, fullfile(graphDir, 'energie.eps'), 'ContentType', 'vector')

% Integral u
fig = figure('Units','centimeters','Position',[5 5 24 12]);
surf(X,Y,integrals_u')
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('DoF Mesh')
ylabel('DoF Circle')
zlabel('Integral u')
shading interp
colorbar
exportgraphics(fig, fullfile(graphDir, 'integral_u_surface.eps'), 'ContentType', 'vector')

% Integral Normalenableitung
fig = figure('Units','centimeters','Position',[5 5 24 12]);
surf(X,Y,integrals_dudn')
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('DoF Mesh')
ylabel('DoF Circle')
zlabel('Integral du/dn')
shading interp
colorbar
exportgraphics(fig, fullfile(graphDir, 'integral_dudn_surface.eps'), 'ContentType', 'vector')

% Laufzeit Plot absolut
idx_circle_h = 5;
gesamtlaufzeit = sw_meshing + sw_solvePDE + sw_solveIntegrals(:, idx_circle_h);

fig = figure('Units','centimeters','Position',[5 5 24 12]);
hold on
    semilogy(num_of_Nodes, sw_meshing + sw_solvePDE + sw_solveIntegrals(:, idx_circle_h))
    semilogy(num_of_Nodes, sw_meshing)
    semilogy(num_of_Nodes, sw_solvePDE)
    semilogy(num_of_Nodes, sw_solveIntegrals(:, idx_circle_h))
hold off
xlabel('DoF Mesh')
ylabel('time [s]')
legend('Gesamtlaufzeit', 'Meshing', 'Lösen der PDE', 'Integrale')
set(gca,'XScale','log')
set(gca,'YScale','log')
exportgraphics(fig, fullfile(graphDir, 'laufzeit_absolut.eps'), 'ContentType', 'vector')

% Laufzeit Plot relativ
fig = figure('Units','centimeters','Position',[5 5 24 12]);
hold on
    semilogy(num_of_Nodes, sw_meshing ./ gesamtlaufzeit)
    semilogy(num_of_Nodes, sw_solvePDE  ./ gesamtlaufzeit)
    semilogy(num_of_Nodes, sw_solveIntegrals(:, idx_circle_h) ./ gesamtlaufzeit)
hold off
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('DoF Mesh')
ylabel('Anteil an der Gesamtlaufzeit')
legend('Meshing','Lösen der PDE', 'Integrale')
exportgraphics(fig, fullfile(graphDir, 'laufzeit_anteile.eps'), 'ContentType', 'vector')

% Laufzeit Integrale Surface
fig = figure('Units','centimeters','Position',[5 5 24 12]);
surf(X,Y,sw_solveIntegrals')
xlabel('DoF Mesh')
ylabel('DoF Circle')
zlabel('time [s]')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'ZScale','log')
shading interp
colorbar
exportgraphics(fig, fullfile(graphDir, 'laufzeit_integrale_surface.eps'), 'ContentType', 'vector')

% Konvergenz I u und I du/dn über DoF Mesh
fig = figure('Units','centimeters','Position',[5 5 24 12]);

subplot(2,1,1)
plot(num_of_Nodes, integrals_u(:,end), 'LineWidth', 1.5)
xticklabels('')
ylabel('Integral u')
set(gca,'XScale','log')

subplot(2,1,2)
plot(num_of_Nodes, integrals_dudn(:,end), 'LineWidth', 1.5)
xlabel('DoF Mesh')
ylabel('Integral du/dn')
set(gca,'XScale','log')

exportgraphics(fig, fullfile(graphDir, 'konvergenz_integrale_mesh.eps'), 'ContentType', 'vector')

% Konvergenz I u und I du/dn über DoF Circle
fig = figure('Units','centimeters','Position',[5 5 24 12]);

subplot(2,1,1)
plot(num_of_C_Nodes, integrals_u(end,:), 'LineWidth', 1.5)
xticklabels('')
ylabel('Integral u')
set(gca,'XScale','log')

subplot(2,1,2)
plot(num_of_C_Nodes, integrals_dudn(end,:), 'LineWidth', 1.5)
xlabel('DoF Circle')
ylabel('Integral du/dn')
set(gca,'XScale','log')

exportgraphics(fig, fullfile(graphDir, 'konvergenz_integrale_circle.eps'), 'ContentType', 'vector')

% # Knoten je hmax
fig = figure('Units','centimeters','Position',[5 5 24 12]);
plot(mesh_sizes, num_of_Nodes, 'LineWidth', 1.5)
xlabel('h mesh')
ylabel('# Knoten')
set(gca,'YScale','log')
exportgraphics(fig, fullfile(graphDir, 'knoten_je_hmesh.eps'), 'ContentType', 'vector')


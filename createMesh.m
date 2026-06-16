function [coordinates, elements, boundary, dirichlet, neumann, hMesh] = createMesh(plotMesh)

if nargin < 1
    plotMesh = false;
end

%% Punkte A, L, K, J, I, H, G, F, E, D, D', C, B

coordinates0 = [
    0,  0;    % A
    1,  5;    % L
    5,  5;    % K
    5,  2;    % J
    9,  0;    % I
    9,  3;    % H
   10,  3;    % G
   10,  0;    % F
   11,  0;    % E
   11, -3;    % D
   10, -3;    % D'
    7, -3;    % C
    5, -3;    % B
];

r_AB = sqrt(34)/2;

g = [
%  A→L  L→K  K→J  J→I  I→H  H→G  G→F  F→E  E→D  D→D' C→D' C→B  A→B   KreisM    KreisN
    2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   1,   2,   1,   1,   1,   1,   1;
    0,   1,   5,   5,   9,   9,  10,  10,  11,  11,   7,   7,   0, 2.5, 2.5, 8.5, 8.5;
    1,   5,   5,   9,   9,  10,  10,  11,  11,  10,  10,   5,   5, 2.5, 2.5, 8.5, 8.5;
    0,   5,   5,   2,   0,   3,   3,   0,   0,  -3,  -3,  -3,   0,-0.5,-2.5,-2.5,-3.5;
    5,   5,   2,   0,   3,   3,   0,   0,  -3,  -3,  -3,  -3,  -3,-2.5,-0.5,-3.5,-2.5;
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   1,   0,   0,   0,   0;
    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   1,   0,   1,   1,   1,   1;
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 8.5,   0, 2.5, 2.5, 2.5, 8.5, 8.5;
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   0,-1.5,-1.5,-1.5,  -3,  -3;
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 1.5,   0,r_AB,   1,   1, 0.5, 0.5;
];

%% Mesh erzeugen

hMesh = 0.3;

[p,e,t] = initmesh(g, 'hmax', hMesh);

elements = t(1:3,:)';
coordinates = p';

%% Randsegmente definieren

boundarys = {
    [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], ...
    [11], [12], [13], [14, 15], [16, 17]
};

boundary = cell(size(boundarys));

for k = 1:length(boundarys)

    idx = find(ismember(e(5,:), boundarys{k}));
    boundary{k} = e(1:2,idx)';

    idxReverse = find(e(6,idx) == 0);
    boundary{k}(idxReverse,[2,1]) = boundary{k}(idxReverse,:);

end

%% Dirichlet- und Neumann-Rand

dirichlet = [
    boundary{6};
    boundary{13};
    boundary{14};
    boundary{15}
];

neumann = [
    boundary{2};
    boundary{3};
    boundary{4};
    boundary{5};
    boundary{7};
    boundary{8};
    boundary{9};
    boundary{10};
    boundary{11};
    boundary{12}
];

%% Optionaler Mesh-Plot

if plotMesh

    figure;

    trisurf(elements, coordinates(:,1), coordinates(:,2), ...
        0*coordinates(:,2), 'edgecolor', 'k', 'facecolor', 'none');

    view(2)
    axis equal
    axis off

    dirichlet_idx = [6, 13, 14, 15];
    neumann_idx = setdiff(1:length(boundary), dirichlet_idx);

    hold on

    h1 = [];
    h2 = [];

    for k = dirichlet_idx
        h1 = plot(reshape(p(1, boundary{k}),[],2)', ...
                  reshape(p(2, boundary{k}),[],2)', ...
                  'r', 'linewidth', 3);
    end

    for k = neumann_idx
        h2 = plot(reshape(p(1, boundary{k}),[],2)', ...
                  reshape(p(2, boundary{k}),[],2)', ...
                  'b', 'linewidth', 3);
    end

    hold off

    view(2)
    axis equal
    axis off

    legend([h1(1), h2(1)], ...
        {'Dirichlet Rand', 'Neumann Rand'}, ...
        'Location', 'best');

end

end
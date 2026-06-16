function [x,energy] = solveLaplace(coordinates,elements,dirichlet,neumann,f,g,uD)

nE = size(elements,1);
nC = size(coordinates,1);
x = zeros(nC,1);

%*** First vertex of elements and corresponding edge vectors
c1  = coordinates(elements(:,1),:);
d21 = coordinates(elements(:,2),:) - c1;
d31 = coordinates(elements(:,3),:) - c1;

%*** Vector of element areas 4*|T|
area4 = 2*(d21(:,1).*d31(:,2) - d21(:,2).*d31(:,1));

%*** Assembly of stiffness matrix
I = reshape(elements(:,[1 2 3 1 2 3 1 2 3])',9*nE,1);
J = reshape(elements(:,[1 1 1 2 2 2 3 3 3])',9*nE,1);

a = (sum(d21.*d31,2)./area4)';
b = (sum(d31.*d31,2)./area4)';
c = (sum(d21.*d21,2)./area4)';

A = [-2*a+b+c; a-b; a-c; a-b; b; -a; a-c; -a; c];
A = sparse(I,J,A(:));

%*** Prescribe values at Dirichlet nodes
dirichlet = unique(dirichlet);
x(dirichlet) = feval(uD,coordinates(dirichlet,:));

%*** Assembly of right-hand side
fsT = feval(f, c1 + (d21 + d31)/3);
b = accumarray(elements(:), repmat(area4.*fsT/12,3,1), [nC 1]) - A*x;

if ~isempty(neumann)
    cn1 = coordinates(neumann(:,1),:);
    cn2 = coordinates(neumann(:,2),:);
    gmE = feval(g, (cn1 + cn2)/2);
    
    edgeLength = sqrt(sum((cn2 - cn1).^2,2));
    b = b + accumarray(neumann(:), ...
        repmat(edgeLength.*gmE/2,2,1), [nC 1]);
end

%*** Computation of P1-FEM approximation
freenodes = setdiff(1:nC, dirichlet);
x(freenodes) = A(freenodes,freenodes) \ b(freenodes);

%*** Compute energy ||grad(uh)||^2 of discrete solution
energy = x' * A * x;

end
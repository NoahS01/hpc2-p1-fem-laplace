function coeffElem = computeElementCoefficients(coordinates, elements, u)

numElem = size(elements,1);
coeffElem = zeros(numElem,3);

for k = 1:numElem

    nodes = elements(k,:);
    P = coordinates(nodes,:);
    U = u(nodes);

    % u_h(x,y) = a + b*x + c*y
    M = [ones(3,1), P(:,1), P(:,2)];

    coeff = M \ U;

    coeffElem(k,:) = coeff';

end

end
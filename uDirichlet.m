function W = uDirichlet(x, coords5, coords3, coords7)

tol = 1e-10;

n = size(x,1);
W = zeros(n,1);

for i = 1:n

    if any(vecnorm(coords5 - x(i,:), 2, 2) < tol)
        W(i) = 5;

    elseif any(vecnorm(coords3 - x(i,:), 2, 2) < tol)
        W(i) = 3;

    elseif any(vecnorm(coords7 - x(i,:), 2, 2) < tol)
        W(i) = 7;
    end

end

end
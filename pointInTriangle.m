function inside = pointInTriangle(p, P)

A = [
    1      1      1;
    P(1,1) P(2,1) P(3,1);
    P(1,2) P(2,2) P(3,2)
];

rhs = [1; p(1); p(2)];

lambda = A \ rhs;

tol = 1e-12;

inside = all(lambda >= -tol) && all(lambda <= 1 + tol);

end
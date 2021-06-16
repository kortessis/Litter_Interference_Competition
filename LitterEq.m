function Leq = LitterEq(lambda, alpha, beta, bdratio, BTdratio)


if beta == 0
    Leq = bdratio*(lambda - 1)/alpha + BTdratio;
elseif alpha/bdratio == 0
    Leq = (lambda-1)/beta;
else
    X = -(1/alpha)*bdratio - 1/beta + BTdratio;
    Y = 1/(beta)*(bdratio*(lambda-1)/alpha + BTdratio);
    Leq = 0.5*X + sqrt(0.25*X.^2 + Y);
end
end

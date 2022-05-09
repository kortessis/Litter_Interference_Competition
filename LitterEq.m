function Leq = LitterEq(R0, alpha, beta, bdratio, BTdratio)


if beta == 0
    Leq = bdratio*(R0 - 1)./alpha + BTdratio;
elseif alpha/bdratio == 0
    Leq = (R0-1)/beta;
else
    X = -(1./alpha)*bdratio - 1/beta + BTdratio;
    Y = 1/(beta)*(bdratio*(R0-1)./alpha + BTdratio);
    Leq = 0.5*X + sqrt(0.25*X.^2 + Y);
end
end

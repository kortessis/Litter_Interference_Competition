function Leq = LitterEq(R0, alpha, beta, bdratio, BTdratio)

% This function evaluates the expression for equilibrium litter density
% found in tables S1 and S2 of the supplement. 

% It takes 5 values as inputs and outputs the equilibirum litter value.
% The 5 inputs are
% 
%     1. R0 <- This is the reproductive number for the species, given by 
%     equations (7) and (8) in the main text. 

%     2. alpha <- This is the effective per-capita competitive effect of a 
%     species. For the annual, input \alpha_A. 
%     **********For the perennial, input
%            \alpha_P*[1 + \gamma*(1-p_P)/p_s],
%      as this is the effective competiitve effect of the perennial. 

%     3. beta <- This is the sensitivity of establishment to litter, i.e.,
%     the \beta's in the expression for E (eq. 6 of the main text). 

%     4. bdratio <- This is the ratio of litter production per-individual
%     to the decomposition fraction. For the annual, it is simply, 
%           b_A/d
%     ***********for the perennial, input
%           (b_P/d)*[\delta + (1-p_P)/p_S].
%
%     5. BTdratio = b_T/d is the ratio of litterfall from trees divided by
%     the decomposition fraction. 
%
if beta == 0
    Leq = bdratio*(R0 - 1)/alpha + BTdratio;
elseif alpha/bdratio == 0
    Leq = (R0-1)/beta;
else
    X = -(1/alpha)*bdratio - 1/beta + BTdratio;
    Y = 1/(beta)*(bdratio*(R0-1)/alpha + BTdratio);
    Leq = 0.5*X + sqrt(0.25*X.^2 + Y);
end
end

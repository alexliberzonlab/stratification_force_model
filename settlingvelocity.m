function ws = settlingvelocity(rhop,rhof,g,d,nu)

Cd = @(Re)  0.25 + (24./Re) + (6./(1+Re.^(0.5)));  % drag law

Ap=(pi*d^2)/4;              % surface of the sphere - m2
Vp=(pi*d^3)/6;              % volume  of the sphere - m3

ws = fzero( @(wp) (rhop-rhof) * Vp * g ...
                   - 0.5 * Cd(wp*d/nu)* rhof * Ap * abs(wp) * wp, ...
           [1e-11 1]);
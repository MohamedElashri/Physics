function [spin,Mag] = Intialization(L, IC)

if IC == 1
    spin = ones(L, L);
elseif IC == 2
    spin = -ones (L, L);
elseif IC == 3
    spin = sign(0.5 - rand(L, L));
end


global Mag;
global En;
global EnMean;
global MagMean;

En = 0;
Mag = 0;

for x=1:L
    for y =1 : L
        N = Neighbor(L, x, y);
        Mag = Mag + spin(x,y);    % -------------------------initial mag
        En = En - spin(x,y)*(spin(N(4),y)+spin(x,N(2)));% --initial energy
    end
end
En;
Mag;
EnMean  = En / L^2;
MagMean = Mag / L^2;
end
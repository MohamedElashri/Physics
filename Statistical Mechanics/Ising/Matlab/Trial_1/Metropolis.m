function [spin] = Metropolis(spin, T, En, Mag, L, n, EnMean, MagMean)
 global En;
 global Mag;
 global EnMean;
 global MagMean;
 global En2Mean;
 global Mag2Mean;
 global Check;
 j = 1;
 E = 0; M = 0; En2 = 0; Mag2 = 0;
 y = 1;
for i = 1:n
    linearIndex = randi(numel(spin));
    [row, col]  = ind2sub(size(spin), linearIndex);
    N = Neighbor(L, row, col);
    dE = 2 * spin(row, col) *(spin(N(1),col)+ spin(row,N(2))+ spin(row,N(3)) + spin(N(4),col)); 
    
    if dE >= 0
        prob = exp(-dE /T);
        x = rand();
         if  x<= prob
             spin(row, col) = - spin(row, col);
             En = En + dE;
             Mag = Mag + 2 * spin(row,col);
         end
    else
         spin(row, col) = - spin(row, col); 
         En = En + dE;                      %--sum of energy
         Mag = Mag + 2 * spin(row,col);
    end
   
   if (Check == 1)
       E = E + En/(L^2);
       M = M + Mag/(L^2);
       En2 = En2 + (En/(L^2))^2;        %---(En.Mean/(L^2))^2
       Mag2 = Mag2 + (Mag/(L^2))^2;     %---(Mag.Mean/(L^2))^2
   end
    if (T==2) || (T==2.5) || ((T == 1) && (Check == 0))
    if mod(i,1000) == 0 
        Energy(j) = En/(L^2);               %----energy per spin
        Magnetization(j) = Mag/(L^2);       %----mag per spin
        if ((mod(j,1000)==0) && (j <= 9000))             
        end
        j= j+1;
    end
    end
end
if (Check == 1)
    EnMean = E/n;                        %-----mean energy
    MagMean = M/n;                     %-----mean mag  
    En2Mean = En2/n;
    Mag2Mean = Mag2/n;
    
end
 
end


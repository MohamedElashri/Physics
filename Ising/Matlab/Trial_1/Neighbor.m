function [Neighbors]  = Neighbor(L,x,y)
above = mod(x - 1 - 1, L) + 1;
below = mod(x + 1 - 1, L) + 1;
left  = mod(y - 1 - 1, L) + 1;
right = mod(y + 1 - 1, L) + 1;
Neighbors = [above ,right, left, below];
end

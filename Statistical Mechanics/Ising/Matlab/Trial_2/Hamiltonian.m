function Emean = Hamiltonian(spin, J)

sumOfNeighbors =  circshift(spin, [ 0  1]) +circshift(spin, [ 0 -1]) + circshift(spin, [ 1  0]) + circshift(spin, [-1  0]);
Em = - J *0.5* spin .* sumOfNeighbors;
Emean = sum(Em(:)) / numel(spin);
end

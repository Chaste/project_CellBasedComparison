% Script to calculate the connectivity matrix of a cell population
% loads the adjacency matrix from a file.

clear 
close all

load /tmp/ozzy/testoutput/CellSorting/Potts/results_from_time_0/cellpopulationadjacency.dat

time = adjacencymatrices(:,1);
num_cells = adjacencymatrices(:,2);


for i=1:length(time)
    adjacency_full = reshape(double(0<adjacencymatrices(i,3:end)),num_cells(i),num_cells(i));
    adjacency_links = reshape(double(3==adjacencymatrices(i,3:end)),num_cells(i),num_cells(i));
    adjacency_common = reshape(double((1==adjacencymatrices(i,3:end)) + (2==adjacencymatrices(i,3:end))),num_cells(i),num_cells(i));
    
%     spy(adjacency_full)
%     spy(adjacency_links)
%     spy(adjacency_common)
   
    max_eig_full(i) = max(eig(adjacency_full));
    num_neighbours_full(i) = sum(sum(adjacency_full));
    
    max_eig_links(i) = max(eig(adjacency_links));
    num_neighbours_links(i) = sum(sum(adjacency_links));
    
    max_eig_common(i) = max(eig(adjacency_common));
    num_neighbours_common(i) = sum(sum(adjacency_common));

end
    
figure;
hold on
plot(time,num_neighbours_full,'k');
plot(time,num_neighbours_links,'r');
plot(time,num_neighbours_common,'b');
plot(time,num_neighbours_common+num_neighbours_links,'kd');

title('num neighbours');
legend('full','links','common', 'links+common');

figure
hold on
plot(time,max_eig_full,'k');
plot(time,max_eig_links,'r');
plot(time,max_eig_common,'b');
legend('full','links','common');
title('max eig')


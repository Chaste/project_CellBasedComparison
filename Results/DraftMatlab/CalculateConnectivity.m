% Script to calculate the connectivity matrix of a cell population
% loads the adjacency matrix from a file.
function PlotComparison()
close all
clear

steady_state = 5;

Type = 'NonProlif'
%Type = 'Prolif'

[time_vertex,num_neighbours_vertex,max_eig_vertex] = CalculateAdjacency([Type '/Vertex']);
[time_mesh,num_neighbours_mesh,max_eig_mesh] = CalculateAdjacency([Type '/MeshGhost']);
[time_node,num_neighbours_node,max_eig_node] = CalculateAdjacency([Type '/Node']);
[time_potts,num_neighbours_potts,max_eig_potts] = CalculateAdjacency([Type '/Potts']);

figure;
hold on
plot(time_vertex-steady_state,num_neighbours_vertex(:,2)./num_neighbours_vertex(:,1),'k');
plot(time_mesh-steady_state,num_neighbours_mesh(:,2)./num_neighbours_mesh(:,1),'k--');
plot(time_node-steady_state,num_neighbours_node(:,2)./num_neighbours_node(:,1),'k-.','linewidth',2.0);
plot(time_potts-steady_state,num_neighbours_potts(:,2)./num_neighbours_potts(:,1),'k:','linewidth',2.0);

xlabel('time');
ylabel('Num Neighbour Ratio');
h = legend('Vertex','Mesh','Node','Potts');
set(h, 'Location', 'northeast')
title(['Num Neighbour Ratio: ', Type]);
saveaspngandeps(-1,[Type 'Figs/NumNeighbourRatio'], 7, 7/5, 9)

figure;
hold on
plot(time_vertex-steady_state,max_eig_vertex(:,2)./max_eig_vertex(:,1),'k');
plot(time_mesh-steady_state,max_eig_mesh(:,2)./max_eig_mesh(:,1),'k--');
plot(time_node-steady_state,max_eig_node(:,2)./max_eig_node(:,1),'k-.','linewidth',2.0);
plot(time_potts-steady_state,max_eig_potts(:,2)./max_eig_potts(:,1),'k:','linewidth',2.0);

xlabel('time');
ylabel('Max Eig Ratio');
h = legend('Vertex','Mesh','Node','Potts');
set(h, 'Location', 'northeast')
title(['Max Eig Ratio: ', Type]);
saveaspngandeps(-1,[Type 'Figs/MaxEigRatio'], 7, 7/5, 9)


function [time,num_neighbours, max_eig] = CalculateAdjacency(Directory)

adjacenciesfile = [Directory,'/results_from_time_5/adjacencymatrices.dat']

adjacenciesdata = LoadNonConstantLengthData(adjacenciesfile);
time = zeros(length(adjacenciesdata),1);
num_cells = zeros(length(adjacenciesdata),1);

for i=1:length(adjacenciesdata)
    adjacency = adjacenciesdata{i};
    time(i) = adjacency(1);
    num_cells(i) = adjacency(2);
    
    adjacency_full = reshape(double(0<adjacency(3:end)),num_cells(i),num_cells(i));
    adjacency_links = reshape(double(3==adjacency(3:end)),num_cells(i),num_cells(i));
    %adjacency_common = reshape(double((1==adjacency(3:end)) + (2==adjacency(3:end))),num_cells(i),num_cells(i));
    
%     spy(adjacency_full)
%     spy(adjacency_links)
%     spy(adjacency_common)
   
    max_eig(i,1) = max(eig(adjacency_full));
    num_neighbours(i,1) = sum(sum(adjacency_full));
    
    max_eig(i,2) = max(eig(adjacency_links));
    num_neighbours(i,2) = sum(sum(adjacency_links));
    
%     max_eig(i,3) = max(eig(adjacency_common));
%     num_neighbours(i,3) = sum(sum(adjacency_common));
end
    

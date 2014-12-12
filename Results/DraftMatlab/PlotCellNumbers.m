function PlotCellNumbers()

clear; close all;

[time_vertex, num_cells_vertex] = LoadCellNumbers('Vertex');
[time_mesh, num_cells_mesh] = LoadCellNumbers('Mesh');
[time_node, num_cells_node] = LoadCellNumbers('Node');
[time_potts, num_cells_potts] = LoadCellNumbers('Potts');

fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
axes('Parent',fig,'FontSize',24);
plot(time_vertex,num_cells_vertex(:,1),'k')
hold on
plot(time_mesh,num_cells_mesh(:,1),'k--')
plot(time_node,num_cells_node(:,1),'k-.','linewidth',2.0)
plot(time_potts,num_cells_potts(:,1),'k:','linewidth',2.0)

% hold on
% plot(MeshFractionalLengths.Var1-steady_state,MeshFractionalLengths.Var2,'k')
% plot(VertexFractionalLengths.Var1-steady_state,NodeFractionalLengths.Var2,'k:')
% plot(VertexFractionalLengths.Var1-steady_state,PottsFractionalLengths.Var2/4,'r')

xlabel('time');
ylabel('Num Cells');
h = legend('Vertex','Mesh','Node','Potts');
set(h, 'Location', 'northwest')
title('Number of Cells')
saveaspngandeps(-1,['NumCells'], 7, 7/5, 9)


fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
axes('Parent',fig,'FontSize',24);
plot(time_vertex,num_cells_vertex(:,2),'k')
hold on
plot(time_mesh,num_cells_mesh(:,2),'k--')
plot(time_node,num_cells_node(:,2),'k-.','linewidth',2.0)
plot(time_potts,num_cells_potts(:,2),'k:','linewidth',2.0)

xlabel('time');
ylabel('Num Cells');
h = legend('Vertex','Mesh','Node','Potts');
set(h, 'Location', 'northwest')
title('Number of LabeledCells')
saveaspngandeps(-1,['NumLabeledCells'], 7, 7/5, 9)

function [time, num_cells] = LoadCellNumbers(Directory)

CellTypesFile = [Directory '/results_from_time_0/results.vizcelltypes']

CellTypesData = LoadNonConstantLengthData(CellTypesFile);

for i=1:length(CellTypesData)
    data = CellTypesData{i};
    time(i) = data(1);
    data(1) = [];
    num_cells(i,1) =length(data)
    num_cells(i,2) =sum(data==1)
end

    


% 
% 
% fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
% axes('Parent',fig,'FontSize',24);
% plot(VertexFractionalLengths.Var1-steady_state,VertexFractionalLengths.Var2./VertexFractionalLengths.Var3,'k--')
% hold on
% plot(MeshFractionalLengths.Var1-steady_state,MeshFractionalLengths.Var2./MeshFractionalLengths.Var3,'k')
% plot(VertexFractionalLengths.Var1-steady_state,NodeFractionalLengths.Var2./NodeFractionalLengths.Var3,'k:')
% plot(VertexFractionalLengths.Var1-steady_state,PottsFractionalLengths.Var2./PottsFractionalLengths.Var3,'r')
% 
% 
% xlabel('time');
% ylabel('Fractional Ratio');
% h = legend('Vertex','Mesh','Node','Potts');
% % set(h, 'Location', 'northwest')
% saveaspngandeps(-1,[Type 'Figs/FractionalLengthRatio'], 7, 7/5, 9)
% 
% 

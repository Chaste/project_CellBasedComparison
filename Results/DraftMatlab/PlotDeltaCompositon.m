% Code to load data averaged across a crypt and plot comparisons
function PlotComparison()
close all
clear

Type = 'NonProlif'
%Type = 'Prolif'


[time_vertex,num_vertex,high_delta_vertex,low_delta_vertex] = CalculateRatios([Type '/Vertex']);
%[time_mesh,num_mesh,high_delta_mesh,low_delta_mesh] = CalculateRatios([Type '/Mesh']);
[time_mesh_ghost,num_mesh_ghost,high_delta_mesh_ghost,low_delta_mesh_ghost] = CalculateRatios([Type '/MeshGhost']);
[time_node,num_node,high_delta_node,low_delta_node] = CalculateRatios([Type '/Node']);
[time_potts,num_potts,high_delta_potts,low_delta_potts] = CalculateRatios([Type '/Potts']);

fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
axes('Parent',fig,'FontSize',24);
plot(time_vertex,(high_delta_vertex+low_delta_vertex)./num_vertex,'k')
hold on
%plot(time_mesh,(high_delta_mesh+low_delta_mesh)./num_mesh,'k')
plot(time_mesh_ghost,(high_delta_mesh_ghost+low_delta_mesh_ghost)./num_mesh_ghost,'k--')
plot(time_node,(high_delta_node+low_delta_node)./num_node,'k-.','linewidth',2.0)
plot(time_potts,(high_delta_potts+low_delta_potts)./num_potts,'k:','linewidth',2.0)


xlabel('time');
ylabel('Ratio');
h = legend('Vertex','Mesh','Node','Potts');
title(['Patterned ratio: ', Type]);
set(h, 'Location', 'southeast')
saveaspngandeps(-1,[Type 'PatternRatio'], 7, 7/5, 9)


function [time,num_cells,high_delta,low_delta] = CalculateRatios(Directory)

low_delta_threshold = 0.05;
high_delta_threshold = 0.95;

deltasfile = [Directory,'/results_from_time_0/celldeltanotch.dat']

deltasdata = LoadNonConstantLengthData(deltasfile);

num_times = length(deltasdata);

time = zeros(num_times,1);
num_cells = time;
high_delta = time;
low_delta = time;

for i= 1:num_times
    time(i) = deltasdata{i}(1);
    cell_id = deltasdata{i}(2:7:end-6);
    delta = deltasdata{i}(6:7:end-2);
    assert(length(cell_id) == length(delta));
    num_cells(i)=length(cell_id);
    high_delta(i) = length(find(delta>high_delta_threshold));
    low_delta(i) = length(find(delta<low_delta_threshold));
end


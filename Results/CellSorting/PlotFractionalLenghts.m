clear; close all;


fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
axes('Parent',fig,'FontSize',24);
VertexFractionalLengths = load('Vertex/results_from_time_0/fractionallengths.dat')
plot(VertexFractionalLengths(:,1),VertexFractionalLengths(:,2),'k')
hold on
plot(VertexFractionalLengths(:,1),VertexFractionalLengths(:,3),'k--')

print -dpng VertexfractionalLength

fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
axes('Parent',fig,'FontSize',24);

PottsFractionalLengths = load('Potts/results_from_time_0/fractionallengths.dat')
plot(PottsFractionalLengths(:,1),PottsFractionalLengths(:,2),'k')
hold on
plot(PottsFractionalLengths(:,1),PottsFractionalLengths(:,3),'k--')

print -dpng PottsFractionalLength
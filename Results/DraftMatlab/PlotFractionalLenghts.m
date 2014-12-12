clear; close all;

steady_state = 5;


Type = 'NonProlif'
%Type = 'Prolif'

VertexFractionalLengths = readtable([Type '/Vertex/results_from_time_5/fractionallengths.dat'],'Delimiter','\t','HeaderLines',1,'ReadVariableNames',false);
MeshFractionalLengths = readtable([Type '/MeshGhost/results_from_time_5/fractionallengths.dat'],'Delimiter','\t','HeaderLines',1,'ReadVariableNames',false);
NodeFractionalLengths = readtable([Type '/Node/results_from_time_5/fractionallengths.dat'],'Delimiter','\t','HeaderLines',1,'ReadVariableNames',false);
PottsFractionalLengths = readtable([Type '/Potts/results_from_time_5/fractionallengths.dat'],'Delimiter','\t','HeaderLines',1,'ReadVariableNames',false);

fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
axes('Parent',fig,'FontSize',24);
plot(VertexFractionalLengths.Var1-steady_state,VertexFractionalLengths.Var2,'k')
hold on
plot(MeshFractionalLengths.Var1-steady_state,MeshFractionalLengths.Var2,'k--')
plot(VertexFractionalLengths.Var1-steady_state,NodeFractionalLengths.Var2,'k-.','linewidth',2.0)
plot(VertexFractionalLengths.Var1-steady_state,PottsFractionalLengths.Var2/4,'k:','linewidth',2.0)

xlabel('time');
ylabel('Fractional Length');
h = legend('Vertex','Mesh','Node','Potts');
set(h, 'Location', 'northeast')
title(['Fractional Length: ',Type])
saveaspngandeps(-1,[Type 'Figs/FractionalLength'], 7, 7/5, 9)


fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
axes('Parent',fig,'FontSize',24);
plot(VertexFractionalLengths.Var1-steady_state,VertexFractionalLengths.Var2./VertexFractionalLengths.Var3,'k')
hold on
plot(MeshFractionalLengths.Var1-steady_state,MeshFractionalLengths.Var2./MeshFractionalLengths.Var3,'k--')
plot(VertexFractionalLengths.Var1-steady_state,NodeFractionalLengths.Var2./NodeFractionalLengths.Var3,'k-.','linewidth',2.0)
plot(VertexFractionalLengths.Var1-steady_state,PottsFractionalLengths.Var2./PottsFractionalLengths.Var3,'k:','linewidth',2.0)

xlabel('time');
ylabel('Fractional Ratio');
h = legend('Vertex','Mesh','Node','Potts');
set(h, 'Location', 'east')
title(['Fractional:Total Length: ',Type])
saveaspngandeps(-1,[Type 'Figs/FractionalLengthRatio'], 7, 7/5, 9)



clear; close all;
FractionalLengths = load('Vertex/results_from_time_0/fractionallengths.dat')
plot(FractionalLengths(:,1),FractionalLengths(:,2),'k')
hold on
plot(FractionalLengths(:,1),FractionalLengths(:,3),'k')
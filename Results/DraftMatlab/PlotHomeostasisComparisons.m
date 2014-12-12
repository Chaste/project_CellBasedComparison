% Code to load data averaged across a crypt and plot comparisons
function PlotComparison()
close all
clear

crypt_height=12;

Type = 'CI'
%Type='Normal'

[height_vertex, averaged_vertex_velocities,averaged_vertex_areas,averaged_vertex_num_cells, vertex_time, vertex_numbers, vertex_areas, vertex_divisions] = CalculateAverages([Type '/Vertex']);
[height_mesh, averaged_mesh_velocities,averaged_mesh_areas,averaged_mesh_num_cells, mesh_time,mesh_numbers, mesh_areas, mesh_divisions] = CalculateAverages([Type '/Mesh']);
[height_node, averaged_node_velocities,averaged_node_areas,averaged_node_num_cells, node_time,node_numbers, node_areas, node_divisions] = CalculateAverages([Type '/Node']);
[height_potts, averaged_potts_velocities,averaged_potts_areas,averaged_potts_num_cells, potts_time,potts_numbers, potts_areas, potts_divisions] = CalculateAverages([Type '/Potts']);
[height_ca, averaged_ca_velocities,averaged_ca_areas,averaged_ca_num_cells, ca_time,ca_numbers, ca_areas, ca_divisions] = CalculateAverages([Type '/Ca']);
% Plot Averaged Velocities

fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
axes('Parent',fig,'FontSize',24);
plot(height_vertex, averaged_vertex_velocities*10, 'k')
hold on
plot(height_mesh, averaged_mesh_velocities*10, 'k--')
plot(height_node, averaged_node_velocities*10, 'k-.','linewidth',2.0)
plot(height_potts/4, averaged_potts_velocities*10/4, 'k:','linewidth',2.0)%/plot(height_ca, averaged_ca_velocities*10, 'b')
%plot([crypt_height/2 crypt_height/2], [-1.0 11.0], 'k:','linewidth',2.0
axis([0,crypt_height,0,4.0]);
xlabel('y');
ylabel('v_y (\mumhr^{-1})');
h = legend('Vertex', 'Mesh', 'Node','Potts');
set(h, 'Location', 'southeast')
title(['Averaged Velocity: ', Type]);
saveaspngandeps(-1,['Figs' Type '/HomeostasisCryptVelocityComparison'], 7, 7/5, 9)


% Plot averaged cell numbers
fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
axes('Parent',fig,'FontSize',24);
plot(height_vertex, averaged_vertex_num_cells, 'k')
hold on
plot(height_mesh, averaged_mesh_num_cells, 'k--')
plot(height_node, averaged_node_num_cells, 'k-.','linewidth',2.0)
plot(height_potts/4, averaged_potts_num_cells, 'k:','linewidth',2.0)
%plot(height_ca, averaged_ca_num_cells, 'b')

xlabel('y');
ylabel('Proportion of total cell number');
axis([0,crypt_height,0.0,0.2]);
h = legend('Vertex', 'Mesh', 'Node','Potts');
set(h, 'Location', 'northeast')
title(['Averaged Cell Numbers: ', Type]);
saveaspngandeps(-1,['Figs' Type '/HomeostasisCryptAveragedNumberComparison'], 7, 7/5, 9)


% Plot averaged cell area
fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
axes('Parent',fig,'FontSize',24);
plot(height_vertex, averaged_vertex_areas, 'k')
hold on
plot(height_mesh, averaged_mesh_areas, 'k--')
plot(height_node, averaged_node_areas, 'k-.','linewidth',2.0)
plot(height_potts/4, averaged_potts_areas/16, 'k:','linewidth',2.0)
%plot(height_ca, averaged_ca_areas, 'b')

xlabel('y');
ylabel('Averaged cell area');
axis([0,crypt_height,0.6,1.0]);
h = legend('Vertex', 'Mesh', 'Node','Potts');
set(h, 'Location', 'southeast')
title(['Averaged Cell Area: ', Type]);
saveaspngandeps(-1,['Figs' Type '/HomeostasisCryptAveragedAreaComparison'], 7, 7/5, 9)



% Plot Cell numbers
steady_state_time = 0;

fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
axes('Parent',fig,'FontSize',24);
plot(vertex_time(steady_state_time+1:end)-steady_state_time, vertex_numbers(steady_state_time+1:end,4),'k')
hold on
plot(mesh_time(steady_state_time+1:end)-steady_state_time, mesh_numbers(steady_state_time+1:end,4), 'k--')
plot(node_time(steady_state_time+1:end)-steady_state_time, node_numbers(steady_state_time+1:end,4), 'k-.','linewidth',2.0);
plot(potts_time(steady_state_time+1:end)-steady_state_time, potts_numbers(steady_state_time+1:end,4), 'k:','linewidth',2.0)
%plot(ca_time(steady_state_time+1:end)-steady_state_time, ca_numbers(steady_state_time+1:end,4), 'b')

plot(vertex_time(steady_state_time+1:end)-steady_state_time, vertex_numbers(steady_state_time+1:end,3),'k')
plot(mesh_time(steady_state_time+1:end)-steady_state_time, mesh_numbers(steady_state_time+1:end,3), 'k--')
plot(node_time(steady_state_time+1:end)-steady_state_time, node_numbers(steady_state_time+1:end,3), 'k-.','linewidth',2.0);
plot(potts_time(steady_state_time+1:end)-steady_state_time, potts_numbers(steady_state_time+1:end,3), 'k:','linewidth',2.0)
%plot(ca_time(steady_state_time+1:end)-steady_state_time, ca_numbers(steady_state_time+1:end,3), 'b')


xlabel('Time (hours)');
ylabel('Number of cells');
axis([0,1000,0,100]);
h = legend('Vertex', 'Mesh', 'Node','Potts');
set(h, 'Location', 'southeast')
title(['Number of Cells: ', Type]);
saveaspngandeps(-1,['Figs' Type '/HomeostasisCryptNumberComparison'], 7, 7/5, 9)
%title('number');


% Plot Division Dist
fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
bin_width = 2;
bins =0:bin_width:crypt_height;
n_vertex = histc(vertex_divisions,bins); % Get number of elements in each bin
n_mesh = histc(mesh_divisions,bins); % Get number of elements in each bin
n_node = histc(node_divisions,bins); % Get number of elements in each bin
n_potts = histc(potts_divisions,bins); % Get number of elements in each bin
%n_ca = histc(ca_divisions,bins); % Get number of elements in each bin

% rescale so proportion
n_vertex = n_vertex./sum(n_vertex).*bin_width;
n_mesh = n_mesh./sum(n_mesh).*bin_width;
n_node = n_node./sum(n_node).*bin_width;
n_potts = n_potts./sum(n_potts).*bin_width;

plot(bins,n_vertex, 'k')
hold on
plot(bins,n_mesh, 'k--')
plot(bins,n_node, 'k-.','linewidth',2.0)
plot(bins,n_potts, 'k:','linewidth',2.0)
xlabel('y');
ylabel('Proportion of total cell number');
title(['Distribution of Cell Divisions: ', Type]);

h = legend('Vertex', 'Mesh', 'Node','Potts');
set(h, 'Location', 'northeast')

saveaspngandeps(-1,['Figs' Type '/CellDivisionDistributionComparison'], 7, 7/5, 9)


% Plot Area Dist
fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
bin_width = 0.2;
potts_areas = potts_areas/16; %Resale potts
bins =0:bin_width:1;
n_vertex = histc(vertex_areas,bins); % Get number of elements in each bin
n_mesh = histc(mesh_areas,bins); % Get number of elements in each bin
n_node = histc(node_areas,bins); % Get number of elements in each bin
n_potts = histc(potts_areas,bins); % Get number of elements in each bin
%n_ca = histc(ca_areas,bins); % Get number of elements in each bin

bar(bins,[n_vertex;n_mesh;n_node;n_potts]',1.0)
legend('Vertex','Mesh','Node','Potts','Location','NorthWest')
% 
% bar(min(vertex_areas):bin_width:max(vertex_areas),n_vertex./sum(n_vertex)./bin_width,'k');
% hold on 
% bar(min(mesh_areas):bin_width:max(mesh_areas),n_mesh./sum(n_mesh)./bin_width,'r')
% bar(min(node_areas):bin_width:max(node_areas),n_node./sum(n_node)./bin_width,'b')
% bar(min(potts_areas):bin_width:max(potts_areas),n_potts./sum(n_potts)./bin_width,'g')
% %bar(min(ca_areas):bin_width:max(ca_areas),n_ca./sum(n_ca)./bin_width,'y')
xlabel('Area (\mum^2)');
ylabel('Proportion of total cell number');
title(['Distribution of Cell Areas: ', Type]);
saveaspngandeps(-1,['Figs' Type '/CellAreaDistributionComparison'], 7, 7/5, 9)

mean_vertex_area = mean(vertex_areas)
mean_mesh_area = mean(mesh_areas)
mean_node_area = mean(node_areas)
mean_potts_area = mean(potts_areas)
%mean_ca_area = mean(ca_areas)
% % Plotting figures for each of the models individually
% % vertex model
% fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
% axes('Parent',fig,'FontSize',24);
% plot(vertex_time(steady_state_time+1:end)-steady_state_time, vertex_numbers(steady_state_time+1:end,3), 'k')
% hold on
% plot(vertex_time(steady_state_time+1:end)-steady_state_time, vertex_numbers(steady_state_time+1:end,1), 'r')
% plot(vertex_time(steady_state_time+1:end)-steady_state_time, vertex_numbers(steady_state_time+1:end,2), 'b')
% xlabel('Time (hours)');
% ylabel('Number of cells');
% title('Vertex')
% %axis([0,1000,0,400]);
% saveaspngandeps(-1,'Figs/VertexNumberComparison', 7, 7/5, 9)
% %title('number');
% 
% % Mesh model
% fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
% axes('Parent',fig,'FontSize',24);
% plot(mesh_time(steady_state_time+1:end)-steady_state_time, mesh_numbers(steady_state_time+1:end,3), 'k')
% hold on
% plot(mesh_time(steady_state_time+1:end)-steady_state_time, mesh_numbers(steady_state_time+1:end,1), 'r')
% plot(mesh_time(steady_state_time+1:end)-steady_state_time, mesh_numbers(steady_state_time+1:end,2), 'b')
% xlabel('Time (hours)');
% ylabel('Number of cells');
% title('Mesh')
% %axis([0,1000,0,400]);
% saveaspngandeps(-1,'Figs/MeshNumberComparison', 7, 7/5, 9)
% %title('number');
% 
% 
% % Node model
% fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
% axes('Parent',fig,'FontSize',24);
% plot(node_time(steady_state_time+1:end)-steady_state_time, node_numbers(steady_state_time+1:end,3), 'k')
% hold on
% plot(node_time(steady_state_time+1:end)-steady_state_time, node_numbers(steady_state_time+1:end,1), 'r')
% plot(node_time(steady_state_time+1:end)-steady_state_time, node_numbers(steady_state_time+1:end,2), 'b')
% xlabel('Time (hours)');
% ylabel('Number of cells');
% title('Node');
% %axis([0,1000,0,80]);
% saveaspngandeps(-1,'Figs/NodeNumberComparison', 7, 7/5, 9)
% %title('number');
% 
% % Potts model
% fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
% axes('Parent',fig,'FontSize',24);
% plot(potts_time(steady_state_time+1:end)-steady_state_time, potts_numbers(steady_state_time+1:end,3), 'k')
% hold on
% plot(potts_time(steady_state_time+1:end)-steady_state_time, potts_numbers(steady_state_time+1:end,1), 'r')
% plot(potts_time(steady_state_time+1:end)-steady_state_time, potts_numbers(steady_state_time+1:end,2), 'b')
% xlabel('Time (hours)');
% ylabel('Number of cells');
% title('Potts');
% %axis([0,1000,0,80]);
% saveaspngandeps(-1,'Figs/PottsNumberComparison', 7, 7/5, 9)
% 
% %title('number');% Ca model
% fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
% axes('Parent',fig,'FontSize',24);
% plot(ca_time(steady_state_time+1:end)-steady_state_time, ca_numbers(steady_state_time+1:end,3), 'k')
% hold on
% plot(ca_time(steady_state_time+1:end)-steady_state_time, ca_numbers(steady_state_time+1:end,1), 'r')
% plot(ca_time(steady_state_time+1:end)-steady_state_time, ca_numbers(steady_state_time+1:end,2), 'b')
% xlabel('Time (hours)');
% ylabel('Number of cells');
% title('Ca');
% %axis([0,1000,0,80]);
% saveaspngandeps(-1,'Figs/CaNumberComparison', 7, 7/5, 9)
% %title('number');


function [bin_centres , averaged_velocities,averaged_areas,averaged_num_cells,time,num_cells,total_area, total_divisions] = CalculateAverages(Directory)

velocitiesfile = [Directory,'/results_from_time_0/cellvelocities.dat']

velocitiesdata = LoadNonConstantLengthData(velocitiesfile);

areasfile = [Directory,'/results_from_time_0/cellareas.dat']
areasdata = LoadNonConstantLengthData(areasfile);

cellsfile = [Directory,'/results_from_time_0/celltypes.dat']
cellsdata = LoadNonConstantLengthData(cellsfile);

divisionsfile = [Directory,'/results_from_time_0/divisions.dat']
divisionsdata = load(divisionsfile);

step=1.0;
crypt_height = 12;

if strcmpi(Directory,'Normal/Potts')||strcmpi(Directory,'CI/Potts')
    step = step*4
    crypt_height = crypt_height*4;
end

height=[0:step:crypt_height];
bin_centres = height(1:end-1) + step/2.0
averaged_velocities = zeros(1,length(bin_centres));
averaged_num_cells = zeros(1,length(bin_centres));
averaged_areas = zeros(1,length(bin_centres));
averaged_num_cells_2 = zeros(1,length(bin_centres));

num_times = length(velocitiesdata)

total_area =[];

time = zeros(num_times,1);
num_cells = zeros(num_times,4);
for i= 1:num_times
    time(i) = velocitiesdata{i}(1);
    num_cells(i,:)=[cellsdata{i}(2) cellsdata{i}(3) cellsdata{i}(4) cellsdata{i}(2)+cellsdata{i}(3)+cellsdata{i}(4)];
end

steady_state = 200; 

total_divisions = divisionsdata(divisionsdata(:,1)>steady_state,3);

Current_Locations = [0;0];
for i= steady_state:1:num_times
    
    cell_id = velocitiesdata{i}(2:5:end-4);
    x = velocitiesdata{i}(3:5:end-3);
    y = velocitiesdata{i}(4:5:end-2);
    u = velocitiesdata{i}(5:5:end-1);
    v = velocitiesdata{i}(6:5:end);
    
    Previous_Locations = Current_Locations;
    Current_Locations = [cell_id; y];
    
    for j=1:length(cell_id)
        bin = 1;
        for k=1:length(height)-1
            if (y(j)<=height(k+1) && y(j)>height(k))
                bin = k;
            end
        end
        %assert(y(j)<height(k+1))
        
        averaged_velocities(bin) = averaged_velocities(bin) + v(j);
        averaged_num_cells(bin) = averaged_num_cells(bin)+1;
    end        

%         previous_location_index = find(Previous_Locations(1,:)==cell_id(j));
%         if (previous_location_index)
%             location_height = Current_Locations(2,j);
%             velocity = Current_Locations(2,j) - Previous_Locations(2,previous_location_index);
%             
%             bin = 1;
%             %if (location_height > 0.1)
%             for k=1:length(height)-1
%                 if (location_height<=height(k+1) && location_height>height(k))
%                     bin = k;
%                 end
%             end
%             if (location_height>height(k+1))
%                 %pause()
%                 k
%                 location_height
%             end
%             
%             averaged_velocities(bin) = averaged_velocities(bin) + velocity;
%             averaged_num_cells(bin) = averaged_num_cells(bin)+1;
%             
%             %end
%         end
%     end
     
    % NOw calculate averaged area note files sample at different times so
    % may have different numbers of cells
    
    cell_id = areasdata{i}(2:5:end-4);
    y = areasdata{i}(5:5:end-1);
    areas = areasdata{i}(6:5:end);
    
    %Save a vecotr of all the areas for plotting histograms
    total_area = [total_area areas];
    
    for j=1:length(cell_id)
        bin = 1;
        for k=1:length(height)-1
            if (y(j)<=height(k+1) && y(j)>height(k))
                bin = k;
            end
        end
        %assert(y(j)<height(k+1))
        
        averaged_areas(bin) = averaged_areas(bin) + areas(j);
        averaged_num_cells_2(bin) = averaged_num_cells_2(bin)+1;
        
    end
end

 
% for i= steady_state:1:num_times
%
%     x = velocitiesdata{i}(3:5:end-3);
%     y = velocitiesdata{i}(4:5:end-2);
%     u = velocitiesdata{i}(5:5:end-1);
%     v = velocitiesdata{i}(6:5:end);
%     %area =
%
%     assert(length(y)==length(v))
%
%     for cell_index = 1:length(y);
%
%         %if(y(cell_index)>0.1)
%             bin = 1;
%             for k=1:length(height)-1
%                 if (y(cell_index)<=height(k+1) && y(cell_index)>height(k))
%                     bin = k;
%                 end
%             end
%
%             averaged_velocities(bin) = averaged_velocities(bin) + v(cell_index);
%             averaged_num_cells(bin) = averaged_num_cells(bin)+1;
%         %end
%     end
%
% end

averaged_velocities = averaged_velocities./averaged_num_cells;
averaged_areas = averaged_areas./averaged_num_cells_2;
averaged_num_cells = averaged_num_cells./sum(averaged_num_cells);

% Plot Area Dist
% fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
% total_area = total_area*100; % Dimensionalise
%
% bin_width = 5;
%
% n = histc(total_area,min(total_area):bin_width:max(total_area)); % Get number of elements in each bin
% bar(min(total_area):bin_width:max(total_area),n./sum(n)./bin_width,'k');
% xlabel('Area (\mum^2)');
% ylabel('Proportion of total cell number');
% title('(A)');
% hold on
% plot([86.6 86.6], [0.0 0.3./bin_width], 'k:','linewidth',2.0)
% axis([0 160 0 0.3./bin_width]);
%

% %Calculate fits
% pd = fitdist(total_area','gamma')
% x_values = 0:1:160;
% pdf1 = pdf(pd,x_values);
% hold on
% plot(x_values,pdf1,'r','LineWidth',2)
%
% pd = fitdist(total_area','normal')
% x_values = 0:1:160;
% pdf2 = pdf(pd,x_values);
% hold on
% plot(x_values,pdf2,'b','LineWidth',2)
% legend('Data', 'Equilibrium', 'Gamma fit', 'Normal fit')
%
% saveaspngandeps(-1,'Figs/HomeostasisCylindricalAreaDist', 7, 7/5, 8)
%


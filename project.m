clear all; close all; clc;
name = 'AAAAA AAAAA';
id = 'A00000000';
hw = 'project';
format short; format compact;

%% Set up parameters
global r A m rho g Cd Cm dt goal field
r = 0.11; A = pi*r^2; m = 0.4; rho = 1.2; 
g = 9.81; Cd = 0.3; Cm = 0.6; dt = 1/100;

load('goal.mat')
load('field.mat')
%% Task 1: Perform simlation 1-7
tic
disp('Begin task 1')

% Get total number of kicks:
input_param = importdata('input_parameter.txt',' ',5);
[num_kick, ~] = size(input_param.data);

% Solve for trajectories
for kID = 1:7
    fprintf('Solve for kick # %d \n', kID)
    [X0, Y0, Z0, Umag, theta, phi, omgX, omgY, omgZ] = ...
        read_input('input_parameter.txt', kID); 
    [T{kID}, X{kID}, Y{kID}, Z{kID}, U{kID}, V{kID}, W{kID}] =  ... 
        soccer(X0, Y0, Z0, Umag, theta, phi, omgX, omgY, omgZ);
end

% Create figure 1
disp('Plot trajectories # 1 - 7 in figure 1')
color = 'krbgmcy';
defender_color = 'bgmcr';
figure('unit', 'in', 'position', [1 4 14 5]); 
hold on;

% Plot trajectories
for kID = 1:7
    plot3(X{kID}, Y{kID}, Z{kID}, color(kID), 'LineWidth', 2);
    legend_string{kID} = sprintf('kick # %d', kID);
end

% Plot landing location
max_time = 0;
for kID = 1:7
    plot3(X{kID}(end), Y{kID}(end), Z{kID}(end), ...
          [color(kID) 'o'], 'MarkerFacecolor', color(kID),...
          'MarkerEdgeColor', color(kID), 'MarkerSize', 6);
    max_time = max(max_time, T{kID}(end));
end

% Bring in defenders
for nd = 1:5
    [Dx, Dy, Dz] = defender(nd, max_time);
    surf(Dx, Dy, Dz, 'FaceColor', defender_color(nd), 'EdgeColor','none'); 
end

% Bring in field and goal
plot3(field.X, field.Y, field.Z, 'go', 'MarkerSize',2);
plot3(goal.Xpost, goal.Ypost, goal.Zpost, 'k-', 'LineWidth',3);
plot3(goal.Xnet, goal.Ynet, goal.Znet, 'co', 'MarkerSize',2);

title('Magnus effect due to varying \omega_z');
legend(legend_string);
axis([-45, 45, 0, 65, 0, 10]); 
view(-20.5, 45); 
box on; grid on;
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
set(gca, 'Position', [0.1 0.12 0.85 .7]);
set(gca, 'FontSize', 14);


%% Task 2: Perform simlation 8-13
disp('Begin task 2')

% Solve for trajectories
for kID = 8:13
    fprintf('Solve for kick # %d \n', kID)
    [X0, Y0, Z0, Umag, theta, phi, omgX, omgY, omgZ] = ...
        read_input('input_parameter.txt', kID); 
    [T{kID}, X{kID}, Y{kID}, Z{kID}, U{kID}, V{kID}, W{kID}] =  ... 
        soccer(X0, Y0, Z0, Umag, theta, phi, omgX, omgY, omgZ);
end


%% Task 2A: create figure 1
disp('Plot trajectories # 8-13 in figure 2')

figure('unit', 'in', 'position', [1 4 14 5]); 
hold on;
legend_string = [];
% Plot trajectories
for kID = 8:13
    plot3(X{kID}, Y{kID}, Z{kID}, color(kID-7), 'LineWidth', 2);
    legend_string{kID-7} = sprintf('kick # %d', kID);
end

% Plot landing location
for kID = 8:13
    plot3(X{kID}(end), Y{kID}(end), Z{kID}(end), ...
          [color(kID-7) 'o'], 'MarkerFacecolor', color(kID-7),...
          'MarkerEdgeColor', color(kID-7), 'MarkerSize', 6);
end

% Bring in defenders
for nd = 1:5
    [Dx, Dy, Dz] = defender(nd, T{9}(end));
    surf(Dx, Dy, Dz, 'FaceColor', defender_color(nd), 'EdgeColor','none'); 
end

% Bring in field and goal
plot3(field.X, field.Y, field.Z, 'go', 'MarkerSize',2);
plot3(goal.Xpost, goal.Ypost, goal.Zpost, 'k-', 'LineWidth',3);
plot3(goal.Xnet, goal.Ynet, goal.Znet, 'co', 'MarkerSize',2);

title('Geometry effect on final position');
legend(legend_string);
axis([-45, 45, 0, 65, 0, 10]); 
view(-20.5, 45); 
box on; grid on;
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
set(gca, 'Position', [0.1 0.12 0.85 .7]);
set(gca, 'FontSize', 14);




% Task 2B: Plot energy in figure 3:
disp('Plot PE and KE in figure 3')
figure(3); hold on;
for kID = 8:13 
    PE = m*g*Z{kID};
    KE = 0.5*m*(U{kID}.^2 + V{kID}.^2+ W{kID}.^2);
    
    subplot(2,3,kID-7)
    hold on;
    plot(T{kID},KE,'k-','LineWidth',2)
    plot(T{kID},PE,'r--','LineWidth',2)       
    title_string = sprintf('kick # %d',kID);
    title(title_string);
    xlabel('time (s)');
    ylabel('Energy (J)');
    set(gca,'FontSize',14);
    axis([0 2.5 0 350]); box on; grid on;
    if kID == 11, legend('KE','PE'), end
end

    
% Task 2C: Build the data structure for sim_res
disp('Create structure sim_res')
for kID = 8:13
    sim_res(kID-7).kick_ID = kID;
    sim_res(kID-7).final_time = T{kID}(end);
    kmax = find(abs(W{kID}) == min(abs(W{kID})));
    sim_res(kID-7).max_height_location = ...
        [X{kID}(kmax) Y{kID}(kmax) Z{kID}(kmax)];
    sim_res(kID-7).final_location = ...
        [X{kID}(end) Y{kID}(end) Z{kID}(end)];
    sim_res(kID-7).travel_distance = ...
         sum( sqrt(diff(X{kID}).^2 + ...
                   diff(Y{kID}).^2 + ...
                   diff(Z{kID}).^2) ); 
end

% Task 2D: Write data to file
disp('Create report.txt')
fid = fopen('report.txt', 'w');
fprintf(fid,'%s\n%s\n', name, id); 
fprintf(fid,'kick_ID, final_time (s),travel_distance (m)\n');
for n = 1:length(sim_res)
    fprintf(fid,'%2d %15.9e %15.9e\n', ...
            sim_res(n).kick_ID, sim_res(n).final_time, ...
            sim_res(n).travel_distance);
end
fclose(fid);

%% Task 3: Identify the kick with goal
disp('Search for the kicks with goal')
clear legend_string;
figure('unit','in','position',[1 2 10 5]); 
hold on;
cs = 'krbgmcy';
counter_g = 0;
for kID = 14:100
    fprintf('Solve for kick # %d \n',kID)
    [X0, Y0, Z0, Umag, theta, phi, omgX, omgY, omgZ] = ...
        read_input('input_parameter.txt',kID); 
    [Tg,Xg,Yg,Zg,Ug,Vg,Wg] =  ... 
        soccer(X0,Y0,Z0,Umag,theta,phi,omgX,omgY,omgZ);
    if Yg(end) > max(goal.Ypost) 
        counter_g = counter_g + 1;
        fprintf('kick_ID %d is goal %d\n',kID,counter_g);
        kID_goal(counter_g) = kID;
        legend_string{counter_g} = sprintf('kick # %d',kID);
        plot3(Xg,Yg,Zg,cs(counter_g),'LineWidth',2);
        plot3(Xg(end),Yg(end),Zg(end), [cs(counter_g) 'o'], ...
         'MarkerFacecolor',cs(counter_g),'MarkerEdgeColor',cs(counter_g),'MarkerSize',10);
    end
end

for nd = 1:5
    [Dx,Dy,Dz] = defender(nd,Tg(end));
    surf(Dx,Dy,Dz,'FaceColor',defender_color(nd), 'EdgeColor','none'); 
end

plot3(field.X,field.Y,field.Z,'go','MarkerSize',2);
plot3(goal.Xpost,goal.Ypost,goal.Zpost,'k-','LineWidth',3);
plot3(goal.Xnet,goal.Ynet,goal.Znet,'co','MarkerSize',2)
title(['Simulation with a goal: #' sprintf('%d',kID_goal)]);
legend(legend_string);
axis([-10 10 25 65 0 4]); 
view(-20.5,45); box on; grid on;
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
set(gca,'Position',[0.1 0.12 0.85 .7]);
set(gca,'FontSize',14);


%% Assign results to requested variables
sp1a = evalc('help read_input');
sp1b = evalc('help soccer');
sp1c = 'See figure 1';
sp2a = 'See figure 2';
sp2b = 'See figure 3';
sp2c = sim_res(1); sp2d = sim_res(2); sp2e = sim_res(3); 
sp2f = sim_res(4); sp2g = sim_res(5); sp2h = sim_res(6);
sp2i = evalc('type report.txt');
sp3a = 'Enter grade for figure 4';
fprintf('Total runtime: %f s\n',toc);
disp('-----Project completed--------');


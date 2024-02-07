function [T, X, Y, Z, U, V, W] = ...
    soccer(X0, Y0, Z0, Umag0, theta, phi, omgX, omgY, omgZ)
% SOCCER solves for the trajectory of a soccer ballinitially located at 
% (X0, Y0, Z0) with an initial velocity Umag in the directions given by 
% angles theta and phi. The spin rate of the ball is given by omgX, omgY 
% and omgZ. Outputs are time T, position (X,Y,Z) and components of velocity 
% (U,V,W) along the trajectory.
% Call format:
% [T, X, Y, Z, U, V, W] = soccer(X0, Y0, Z0, Umag0, theta, phi, omgX, omgY, omgZ)

global r A m rho g Cd Cm dt goal field

%% Pre-allocation:
T = zeros(1,250);
X = zeros(1,250);
Y = zeros(1,250);
Z = zeros(1,250);
U = zeros(1,250);
V = zeros(1,250);
W = zeros(1,250);

%% Initialize the kick
n = 1;
T(n) = 0;
X(n) = X0;
Y(n) = Y0;
Z(n) = Z0;
U(n) = Umag0*cosd(theta)*sind(phi);
V(n) = Umag0*sind(theta)*sind(phi);
W(n) = Umag0*cosd(phi);

%% Check bound
XFmax = max(field.X) - r;
XFmin = min(field.X) + r;
YFmax = max(field.Y) - r;
YFmin = min(field.Y) + r;
XGmax = max(goal.Xpost);
XGmin = min(goal.Xpost);
ZGmax = max(goal.Zpost);
YNmax = max(goal.Ynet); 

%% Check whether the ball is intially inside the field and ready to shoot
if  XFmin < X(n) && X(n) < XFmax && ...
    YFmin < Y(n) && Y(n) < YFmax && ...
    Z(n) >= r && Umag0 > 0
    advance = true;
    in_goal = false;
else % exiting the function
    advance = false;
    return
end
   
%% Start simulation
while advance
    % Advance governing equations using Euler-Cromer method
    Vmag = sqrt(U(n)^2+V(n)^2+W(n)^2);
    fric_coeff = Cd*rho*A/2/m;
    Magnus_coeff = Cm*rho*A*r/2/m;
    U(n+1) = U(n) - dt * (fric_coeff * Vmag * U(n) ...
                         - Magnus_coeff * (omgY*W(n) - omgZ*V(n)));
    V(n+1) = V(n) - dt * (fric_coeff * Vmag * V(n)...
                         - Magnus_coeff * (omgZ*U(n) - omgX*W(n)));
    W(n+1) = W(n) - dt * (fric_coeff * Vmag * W(n) ...
                         - Magnus_coeff * (omgX*V(n) - omgY*U(n)) + g);
    X(n+1) = X(n) + dt * U(n+1);
    Y(n+1) = Y(n) + dt * V(n+1);
    Z(n+1) = Z(n) + dt * W(n+1);
    T(n+1) = T(n) + dt;
    
    % Compute distance to defenders
    for nd = 1:5
        [Dx, Dy, Dz] = defender(nd,T(n+1));
        dist_2_defender(nd) = ...
            min(min(sqrt((Dx-X(n+1)).^2 + (Dy-Y(n+1)).^2 + (Dz-Z(n+1)).^2)));
    end
    
    % Apply conditions for termination
    if  Z(n+1) <= r                         % hit ground
        break;
    elseif any(dist_2_defender <= r)        % hit defenders
        break;
    elseif X(n+1) <= XFmin                  % hit left bound                
        break;
    elseif X(n+1) >= XFmax                  % hit right bound
        break;
    elseif Y(n+1) <= YFmin                  % hit front bound
        break;
    elseif Y(n+1) >= YFmax                  % hit back bound
        if Z(n+1) >= (ZGmax - r)            % ...+ above goal
            break;
        elseif XFmin <= X(n+1) && X(n+1) <= (XGmin-r)  % ...+ left of goal
            break;
        elseif (XGmax+r) <= X(n+1) && X(n+1) <= XFmax  % ...+ right of goal
            break;
        else % in goal: enforce net condition
           dist_2_net =  min(min(sqrt((goal.Xnet-X(n+1)).^2  ...
               + (goal.Ynet-Y(n+1)).^2 + (goal.Znet-Z(n+1)).^2)));
           if dist_2_net <= r
               break;
           end
        end
    end           
    n = n+1;
end

%% Remove zeros due to pre-allocation:
X(n+1:end) = []; Y(n+1:end) = []; Z(n+1:end) = []; 
U(n+1:end) = []; V(n+1:end) = []; W(n+1:end) = [];
T(n+1:end) = [];
end % function soccer

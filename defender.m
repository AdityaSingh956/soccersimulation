function [Dx, Dy, Dz] = defender(defender_num, time)
% DEFENDER returns the (x,y,z) data used to render a defender at a specific 
% input time. Input defender_num is an integer 1 to 5 denoting the 
% defender number. Defenders 1-4 jump up down in time. 
% Defender 5,the goalkeeper, jumps to catch the ball. 
% Output Dx, Dy, Dz are two-dimensional matrix coordinates 
% used to render the defender surface.
% Call format: [Dx, Dy, Dz] = defender(Defender_num, time)

n = 40; % number of grid points used to generate defender's surface
switch defender_num
    case 1
        % Location of the center of the ellipsoid
        % Translate the ellipsoid up down in time
        XC = 1; YC = 43.5; ZC = 0.9+0.5*sin(0.5*pi/0.8*time)^2;
        % Radii of the ellipsoid
        XR = 0.4; YR = 0.4; ZR = 0.9;
    case 2
        XC = 0; YC = 43.5; ZC = 0.85+0.5*sin(0.5*pi/0.8*time)^2;
        XR = 0.4; YR = 0.4; ZR = 0.85;
    case 3
        XC = -1; YC = 43.5; ZC = 0.85+0.6*sin(0.5*pi/0.8*time)^2;
        XR = 0.4; YR = 0.4; ZR = 0.85;
    case 4
        XC = -2; YC = 43.5; ZC = 0.91+0.75*sin(0.5*pi/0.8*time)^2;
        XR = 0.4; YR = 0.4; ZR = 0.91;
    case 5 % goalkeeper
        XC = 0; YC = 58; ZC = 0;
        XR = 0.4; YR = 0.4; ZR = 0.92; 
    otherwise
        disp('Invalid Defender number');
end
[Dx, Dy, Dz] = ellipsoid(XC,YC,ZC,XR,YR,ZR,n);

% Simulate goalkeeper's jump
if defender_num == 5
    % Rotation ellipsoid by angle alpha
    alpha = 75*time/1.58; 
    Dxn = (Dx)*cosd(alpha) + (Dz)*sind(alpha);
    Dzn = -(Dx)*sind(alpha) + (Dz)*cosd(alpha);

    % Shift ellipsoid to initial location x = -2, z = .92
    % then translate toward positive x and positive z
    Dx = Dxn - 2 + 5*time/1.5;
    Dz = Dzn + 0.92 + 1.1*time/1.5;
end

end %function defender
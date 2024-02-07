function [X0, Y0, Z0, Umag0, theta, phi, omgX, omgY, omgZ] = ...
                                      read_input(input_filename, kick_id)
% READ_INPUT reads user's parameters from file input_filename and output 
% parameters needed to compute the trajectory of the soccer ball.
% Call format: [X0, Y0, Z0, Umag, theta, phi, omgX, omgY, omgZ] = ...
%            read_input(inputfile, kick_id)


param = importdata( input_filename,' ', 5);

if any(kick_id == param.data(:, 1))
    row = find(kick_id == param.data(:, 1));
    X0 = param.data(row, 2);
    Y0 = param.data(row, 3);
    Z0 = param.data(row, 4);
    Umag0 = param.data(row, 5);
    theta = param.data(row, 6);
    phi = param.data(row, 7);
    omgX = param.data(row, 8);
    omgY = param.data(row, 9);
    omgZ = param.data(row, 10);
else
    disp('read_input: kick_id is not valid');
    Xo = NaN;
    Yo = NaN;
    Zo = NaN;
    Umag0 = NaN;
    theta = NaN;
    phi = NaN;
    omgX = NaN;
    omgY = NaN;
    omgZ = NaN;
end

end % function read_input
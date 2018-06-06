function [ Velocity_top_lid ] = CreateLidVelocity( alpha,N )
%CREATEALPHAINPUT creates a parabolic velcoity input on the top lid with
%amplitude alpha

% form the grid if not already done it
persistent Grid 

if isempty(Grid)
    [Grid]=CollocationGrid_q(N);    % the computational grid
end

Velocity_top_lid = alpha*(1-Grid.x.^2).^2;

end


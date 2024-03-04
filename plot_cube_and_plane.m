function plot_cube_and_plane(volumeSize, normalVector, pointOnPlane)
    % volumeSize is a 3-element vector [X Y Z] with the dimensions of the volume
    % normalVector is the normal vector of the plane [a1 a2 a3]
    % pointOnPlane is a point through which the plane passes [x0 y0 z0]
    
    % Define the vertices of the cube
    [X, Y, Z] = deal(volumeSize(1), volumeSize(2), volumeSize(3));
    vertices = [0 0 0; X 0 0; X Y 0; 0 Y 0; 0 0 Z; X 0 Z; X Y Z; 0 Y Z];
    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
    
    % Plot the cube
    figure;
    hold on;
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'none', 'EdgeColor', 'r');
    
    % Define the plane
    % The plane equation is a1*(x-x0) + a2*(y-y0) + a3*(z-z0) = 0
    % We'll create a grid of points and calculate the corresponding Z values
    [xGrid, yGrid] = meshgrid(linspace(0, X, 10), linspace(0, Y, 10));
    zGrid = (-normalVector(1) * (xGrid - pointOnPlane(1)) - ...
             normalVector(2) * (yGrid - pointOnPlane(2))) / normalVector(3) + pointOnPlane(3);
    
    % Plot the plane
    surf(xGrid, yGrid, zGrid, 'FaceColor', 'cyan', 'FaceAlpha', 0.5);
    
    % Set the aspect ratio based on the volume dimensions
    daspect([1 1 1]);
    view(3); % Set to 3D view
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    grid on;
    hold off;
end

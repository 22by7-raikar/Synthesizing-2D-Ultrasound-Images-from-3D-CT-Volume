function plot_clipped_plane_in_volume(niiFile, normalVector, pointOnPlane)

    % Read NIfTI file and get volume size
    nii = load_nii(niiFile);
    volumeData = nii.img;
    volumeSize = size(volumeData);
    [X, Y, Z] = deal(volumeSize(1), volumeSize(2), volumeSize(3));
    
    % Define the vertices and faces of the cube representing the volume
    vertices = [1 1 1; X 1 1; X Y 1; 1 Y 1; 1 1 Z; X 1 Z; X Y Z; 1 Y Z];
    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
    
    % Plot the cube
    figure;
    hold on;
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'none', 'EdgeColor', 'r');
    
    % Define the grid for the plane within the volume
    [xGrid, yGrid] = meshgrid(1:X, 1:Y);
    % Ensure the normal vector is not zero to avoid division by zero
    if normalVector(3) == 0
        normalVector(3) = eps; % Use MATLAB's smallest value to avoid division by zero
    end
    % Calculate zGrid based on the plane equation
    zGrid = (-normalVector(1) * (xGrid - pointOnPlane(1)) - ...
             normalVector(2) * (yGrid - pointOnPlane(2))) / normalVector(3) + pointOnPlane(3);
         
    % Clip the plane to the volume bounds
    insideVolume = (xGrid >= 1) & (xGrid <= X) & ...
                   (yGrid >= 1) & (yGrid <= Y) & ...
                   (zGrid >= 1) & (zGrid <= Z);
    xGrid(~insideVolume) = NaN;
    yGrid(~insideVolume) = NaN;
    zGrid(~insideVolume) = NaN;
    
    % Plot the clipped plane
    surf(xGrid, yGrid, zGrid, 'FaceColor', 'cyan', 'FaceAlpha', 0.5);
    
    % Set the aspect ratio and labels
    daspect([1 1 1]);
    view(3); % Set to 3D view
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    title('Clipped Plane in Volume');
    grid on;
    hold off;
end

function plot_clipped_plane_in_volume(niiFile, normalVector, pointOnPlane)

    % Load the NIfTI file and get the volume data
    nii = load_nii(niiFile);
    volumeData = nii.img;

    % Define the volume size
    [X, Y, Z] = size(volumeData);
    
    % Define the vertices and faces of the cube representing the volume
    vertices = [1 1 1; X 1 1; X Y 1; 1 Y 1; 1 1 Z; X 1 Z; X Y Z; 1 Y Z];
    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
    
    % Plot the cube
    figure;
    hold on;
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'none', 'EdgeColor', 'r');
    
    % Define the grid for the plane
    [xGrid, yGrid] = meshgrid(1:X, 1:Y);
    % Calculate zGrid based on the plane equation
    zGrid = (-normalVector(1) * (xGrid - pointOnPlane(1)) - ...
             normalVector(2) * (yGrid - pointOnPlane(2))) / normalVector(3) + pointOnPlane(3);
    
    % Ensure the calculated Z-grid is within the volume bounds
    validMask = (xGrid >= 1) & (xGrid <= X) & ...
                (yGrid >= 1) & (yGrid <= Y) & ...
                (zGrid >= 1) & (zGrid <= Z);
                
    % Clip the grids to the valid mask
    xGrid = xGrid .* validMask;
    yGrid = yGrid .* validMask;
    zGrid = zGrid .* validMask;
    
    % Use logical indexing to extract the valid slice values
    sliceIndices = sub2ind(size(volumeData), round(zGrid(validMask)), round(yGrid(validMask)), round(xGrid(validMask)));
    sliceValues = volumeData(sliceIndices);
    
    % Initialize the slice matrix with NaNs (to represent the absence of data)
    sliceMatrix = NaN(size(xGrid));
    sliceMatrix(validMask) = sliceValues;
    
    % Plot the slice
    surf(xGrid, yGrid, zGrid, sliceMatrix, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
    colormap(gray); % Use a grayscale colormap
    colorbar; % Show a colorbar
    
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

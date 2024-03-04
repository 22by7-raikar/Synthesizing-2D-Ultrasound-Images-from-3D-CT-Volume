function slice = get_arbitrary_slice(ctdata, normalVector, pointOnPlane,stepSize, voxelSizeX, voxelSizeY, voxelSizeZ)

    % Read NIfTI file
    volumeData = ctdata;
    
    % Extract the normal vector components
    a1 = normalVector(1);
    a2 = normalVector(2);
    a3 = normalVector(3);
    
    % Extract the point on the plane
    x0 = pointOnPlane(1);
    y0 = pointOnPlane(2);
    z0 = pointOnPlane(3);
    
    % Normalize the normal vector
    normNormal = sqrt(a1^2 + a2^2 + a3^2);
    a1 = a1 / normNormal;
    a2 = a2 / normNormal;
    a3 = a3 / normNormal;
    
    % Determine the dimensions of the CT volume
    [X, Y, Z] = size(volumeData);
    
    % Create a grid of y and z coordinates, scaled by their respective voxel sizes
    [Ygrid, Zgrid] = meshgrid((1:Y) * voxelSizeY, (1:Z) * voxelSizeZ);
    
    % Calculate the corresponding x coordinates on the slice, incorporating voxel sizes
    Xgrid = ((a1/a3) * (Zgrid - z0 * voxelSizeZ) + (a1/a2) * (Ygrid - y0 * voxelSizeY)) / voxelSizeX + x0;
    
    % Flatten the grids for interpolation
YgridLinear = Ygrid(:);
ZgridLinear = Zgrid(:);
XgridLinear = Xgrid(:);

% Perform interpolation for all points at once
% Ensure that we're using valid indices
validIdx = (XgridLinear >= 1) & (XgridLinear <= X) & ...
           (YgridLinear >= 1) & (YgridLinear <= Y) & ...
           (ZgridLinear >= 1) & (ZgridLinear <= Z);
sliceValues = interp3(volumeData, XgridLinear(validIdx), YgridLinear(validIdx), ZgridLinear(validIdx), 'nearest');

% Initialize the slice matrix to the correct size
% Ensure that the sizes are integers
sliceSizeY = round(max(Ygrid(:)));
sliceSizeZ = round(max(Zgrid(:)));
slice = zeros(sliceSizeY, sliceSizeZ);

% Place the interpolated values back into the slice matrix
% Reshape Ygrid and Zgrid to match the slice dimensions before using sub2ind
% Make sure Ygrid and Zgrid are also integers
Ygrid = round(Ygrid(validIdx));
Zgrid = round(Zgrid(validIdx));

% Now use sub2ind to convert subscript to index
indices = sub2ind([sliceSizeY, sliceSizeZ], Ygrid, Zgrid);
slice(indices) = sliceValues;
    % Display the slice
    imagesc(slice);
    colormap(gray); % Use a grayscale colormap appropriate for CT images
    axis image; % Set the axis to image scaling
    
    % Set the data aspect ratio based on voxel dimensions
    daspect([1, voxelSizeX/voxelSizeY, voxelSizeZ/voxelSizeX]);
    
    colorbar; % Optionally show a colorbar to indicate the scale
    title('Extracted Slice from 3D CT Volume');
end

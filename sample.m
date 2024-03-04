% Define the point on the plane
pointOnPlane = [50, 50, 50];

% Read NIfTI file
volumeData = niftiread("CTLiver.nii");

% Convert angles to radians
angleX_rad = deg2rad(90);
angleY_rad = deg2rad(90);

% Calculate the normal vector components
a1 = cos(angleX_rad) * cos(angleY_rad); % x-component
a2 = sin(angleX_rad) * cos(angleY_rad); % y-component
a3 = sin(angleY_rad); % z-component

% Normalize the normal vector
normalVector = [a1, a2, a3] / sqrt(a1^2 + a2^2 + a3^2);

% Point on the plane P0 = [x0, y0, z0]
x0 = pointOnPlane(1);
y0 = pointOnPlane(2);
z0 = pointOnPlane(3);

% Determine the dimensions of the CT volume
[X, Y, Z] = size(volumeData);

% Preallocate the slice matrix
slice = zeros(Y, Z);

% Create Y-Z grid
[Y_grid, Z_grid] = meshgrid(1:Y, 1:Z);

% Calculate the X coordinates for the entire grid
X_grid = x0 + (Y_grid - y0) .* (a1 / a2) + (Z_grid - z0) .* (a1 / a3);

% Ensure the X coordinates are within the bounds
X_grid = max(min(X_grid, X), 1);

% Use single interpolation call for the entire slice
slice(:) = interp3(volumeData, X_grid, Y_grid, Z_grid, 'nearest');

% Visualization
imagesc(slice);
colormap('gray');
axis image; % Set the axis to image scaling
colorbar; % Optionally show a colorbar to indicate the scale
title('Extracted Slice from 3D CT Volume');

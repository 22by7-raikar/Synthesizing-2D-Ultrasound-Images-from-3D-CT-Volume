function 

stepSize = 1
    % Read NIfTI file
    volumeData = niftiread('CTLiver.nii');
    
    normalVector = [1,0.866,0.5]
    % Normal vector n = [a1, a2, a3]
    a1 = normalVector(1);
    a2 = normalVector(2);
    a3 = normalVector(3);
    
    pointOnPlane = [512/2,210,241]
    % Point on the plane P0 = [x0, y0, z0]
    x0 = pointOnPlane(1);
    y0 = pointOnPlane(2);
    z0 = pointOnPlane(3);
    
    % Calculate direction cosines
    cosAlpha = a1 / sqrt(a1^2 + a2^2 + a3^2);
    cosBeta = a2 / sqrt(a1^2 + a2^2 + a3^2);
    cosGamma = a3 / sqrt(a1^2 + a2^2 + a3^2);
    
    % Determine the dimensions of the CT volume
    [X, Y, Z] = size(volumeData);
    
    % Create a grid of y and z coordinates
    [Ygrid, Zgrid] = meshgrid(1:Y, 1:Z);
    
    % Calculate the corresponding x coordinates on the slice
    % Using element-wise operations to calculate for entire grid
    Xgrid = x0 + ((Ygrid - y0) / cosBeta) * cosAlpha + ((Zgrid - z0) / cosGamma) * cosAlpha;
    
    % Flatten the grids for interpolation
    Ygrid = Ygrid(:);
    Zgrid = Zgrid(:);
    Xgrid = Xgrid(:);
    
    % Perform interpolation for all points at once
    % Ensure that we're using valid indices
    validIdx = (Xgrid >= 1) & (Xgrid <= X) & (Ygrid >= 1) & (Ygrid <= Y) & (Zgrid >= 1) & (Zgrid <= Z);
    sliceValues = interp3(volumeData, Xgrid(validIdx), Ygrid(validIdx), Zgrid(validIdx), 'nearest');
    
    % Initialize the slice matrix
    slice = zeros(Y, Z);
    
    % Place the interpolated values back into the slice matrix
    slice(sub2ind([Y, Z], Ygrid(validIdx), Zgrid(validIdx))) = sliceValues;
imagesc(slice);
colormap(gray); % Use a grayscale colormap appropriate for CT images
axis image; % Set the axis to image scaling
colorbar; % Optionally show a colorbar to indicate the scale
title('Extracted Slice from 3D CT Volume');

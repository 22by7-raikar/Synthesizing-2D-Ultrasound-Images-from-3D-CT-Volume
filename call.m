% Assuming `ctVolume` is your 3D CT data loaded into MATLAB
angle = 45; % Desired rotation angle in degrees
sliceIndex = 100; % Index of the slice to extract
axis = 'z'; % Rotate around the z-axis for a horizontal slice
ctvolume

% Get the rotated slice
rotatedSlice = getRotatedSlice(ctVolume, angle, sliceIndex, axis);

% Display the rotated slice
imshow(rotatedSlice, []);



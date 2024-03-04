filename = 'CTLiver.nii';
volumeData = niftiread(filename);
normalVector = [1, 1, 1]; % Example normal vector
pointOnPlane = [50, 50, 50]; % Example point on the plane in voxel coordinates
voxelSize = [0.6992, 0.6992, 0.7]; % Voxel dimensions in mm

visualize_slice_with_cube(volumeData, normalVector, pointOnPlane, voxelSize);

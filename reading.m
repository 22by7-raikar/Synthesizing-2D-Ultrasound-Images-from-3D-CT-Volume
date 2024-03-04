function [ctdata,dimensions,voxelSizes,sliceThickness] = reading(filename)
%% Reading the 3D CT Volume and visualizing a randomly selected slice

ctdata = niftiread(filename);
linfo = niftiinfo(filename);

% Extract dimensions
dimensions = linfo.ImageSize;


% Extract voxel sizes
voxelSizes = abs(linfo.PixelDimensions);

numSlices = linfo.ImageSize(3);

% Extract the slice thickness
sliceThickness = abs(linfo.PixelDimensions(3));

% % Display information
% disp('Number of Slices:');
% disp(numSlices);
% 
% disp('Slice Thickness (mm):');
% disp(sliceThickness);
% 
% % Display information
% disp('Dimensions:');
% disp(dimensions);
% 
% disp('Voxel Sizes (mm):');
% disp(voxelSizes);
% 
% % Display unit information
% disp('Units:');
% disp(linfo.SpaceUnits);

end



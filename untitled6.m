niiFile = 'CTLiver.nii';  % Replace with the actual path to your NIfTI file
normalVector = [250, 250, 512];  % Replace with your actual normal vector
pointOnPlane = [50, 50, 50];  % Replace with your actual point on the plane

plot_clipped_plane_in_volume(niiFile, normalVector, pointOnPlane);

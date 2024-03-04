clear
close all
clc

%% Reading the 3D CT volume Data
% Enter the file name of the NiFTi file you want to slice

filename = 'CTLiver.nii';

% Obtain the CT data and unit info of the 3D CT volume
[ctdata,dimensions,voxelSizes,sliceThickness] = reading(filename);

voxelSizeX = voxelSizes(1); % The physical size of a voxel along the X dimension
voxelSizeY = voxelSizes(2); % The physical size of a voxel along the Y dimension
voxelSizeZ = voxelSizes(3); % The physical size of a voxel along the Z dimension

% Volume size
volumeSize = [dimensions(1), dimensions(2), dimensions(3)];

%% Giving the Orientation of the Slice Plane we want to obtain

%Define the normal vector for the desired plane

normalVector = [0.8,1,0.5];

%Give a point on plane you want the normal vector to pass from

pointOnPlane = [512/2,512/2,682/2]; % Set only if you desire a slice at an orientation

depth = 100;    % Set only if the slice to be cut lies in XY/YZ/XZ plane

stepSize = 1;   %Deafult as 1

%% Obtaining the desired slice from 3D CT volume and saving it for post processing

% XY Slice
if normalVector == [0,0,1]
    c = ctdata(:,:,depth);
% XZ Slice
elseif normalVector == [0,1,0]
    c = squeeze(ctdata(depth,:,:));
% YZ Slice
elseif normalVector == [1,0,0]
    c = squeeze(ctdata(:,depth,:));
% Slicing at a desired pose
else
    c = get_arbitrary_slice(ctdata, normalVector, pointOnPlane, stepSize, voxelSizeX, voxelSizeY, voxelSizeZ);
end

%% Save the slice image to a file for post processing
sliceFileName = 'sliceCTImage.jpg';
imwrite(c,sliceFileName);
imshow("sliceCTImage.jpg")


%% Visualization of the Slice in 3D plane 
%Comment them if you donot want to visualize the slice

% Plot the cube
[X, Y, Z] = deal(volumeSize(1), volumeSize(2), volumeSize(3));
vertices = [0 0 0; X 0 0; X Y 0; 0 Y 0; 0 0 Z; X 0 Z; X Y Z; 0 Y Z];
faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];

figure;
hold on;
patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'none', 'EdgeColor', 'r');

% Define the plane
[xGrid, yGrid] = meshgrid(linspace(1, X, size(c, 2)), linspace(1, Y, size(c, 1)));
zGrid = (-normalVector(1) * (xGrid - pointOnPlane(1)) - ...
         normalVector(2) * (yGrid - pointOnPlane(2))) / normalVector(3) + pointOnPlane(3);

% Read the saved slice image
sliceImg = imread(sliceFileName);
sliceImg = imresize(sliceImg, [size(xGrid, 1), size(xGrid, 2)]);
% Superimpose the image onto the plane
h = surf(xGrid, yGrid, zGrid, 'CData', sliceImg, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
alpha(h, 0.5); % Set the transparency to 50%

% Set the aspect ratio based on the volume dimensions
daspect([1 1 1]);
view(3); % Set to 3D view
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('CT Volume with Superimposed Slice');
grid on;
hold off;

% Display the extracted slice for verification
figure;
imagesc(c);
colormap(gray); % Use a grayscale colormap appropriate for CT images
daspect([1, voxelSizeX/voxelSizeY, voxelSizeZ/voxelSizeX]);
colorbar; % Show a colorbar
title('Extracted Slice from 3D CT Volume');

 %% MUST PIPELINE TO SYNTHESIZE US IMAGES FROM CT SLICE
 
%% Select a transducer with GETPARAM

%The structure param contains the tranducer properties.
param = getparam('P4-2v'); % a 2.7-MHz 64-element phased array.

%Add transmit apodization
param.TXapodization = cos(linspace(-1,1,64)*pi/2);

%% Set the transmit delays with TXDELAY

%The left ventricle will be insonified with seven diverging waves 60 degrees wide and tilted at -20 to +20 degrees.

tilt = deg2rad(linspace(-20,20,7)); % tilt angles in rad
txdel = cell(7,1); % this cell will contain the transmit delays     

%Use TXDELAY to calculate the transmit delays for the 7 diverging waves.
for k = 1:7
    txdel{k} = txdelay(param,tilt(k),deg2rad(60));
end

%These are the transmit delays to obtain a 60 degrees wide circular waves steered at -20 degrees
stem(txdel{1}*1e6)
xlabel('Element number')
ylabel('Delays (\mus)')
title('TX delays for a 60{\circ}-wide -20{\circ}-tilted wave')
axis tight square

%% Simulate an acoustic pressure field with PFIELD
% Check what the sound pressure fields look like.

% Define a 100 x 100 polar grid using IMPOLGRID:
[xi,zi] = impolgrid([100 100],15e-2,deg2rad(120),param);

% Simulate an RMS pressure field...
option.WaitBar = false;
P = pfield(xi,0*xi,zi,txdel{1},param,option);

%% ..and display the result:

pcolor(xi*1e2,zi*1e2,20*log10(P/max(P(:))))
shading interp
xlabel('x (cm)')
ylabel('z (cm)')
title('RMS pressure field for a 60{\circ}-wide -20{\circ}-tilted wave')
axis equal ij tight
clim([-20 0]) % dynamic range = [-20,0] dB
cb = colorbar;
cb.YTickLabel{end} = '0 dB';
colormap(hot)

%% Simulate RF signals with SIMUS

% We will now simulate seven series of RF signals. Each series will contain 
% 64 RF signals, as the simulated phased array contains 64 elements. We 
% first create scatterers with GENSCAT from a clipart image of the ventricle stored in heart.jpg

% x, y, and z contain the scatterers' positions. RC contains the reflection coefficients.
% I = rgb2gray(imread('liver100.jpg'));

I = im2gray(imread(sliceFileName));

% Pseudorandom distribution of scatterers (depth is 15 cm)
[x,y,z,RC] = genscat([NaN 15e-2],1540/param.fc,I);

% Take a look at the scatterers. The backscatter coefficients are gamma-compressed to increase the contrast.

scatter(x*1e2,z*1e2,2,abs(RC).^.25,'filled')
colormap([1-hot;hot])
axis equal ij tight
set(gca,'XColor','none','box','off')
title('Scatterers for an abdominal CT slice')


% Simulate the seven series of RF signals with SIMUS. 
% The RF signals will be sampled at 4 $\times$ center frequency.

RF = cell(7,1); % this cell will contain the RF series
param.fs = 4*param.fc; % sampling frequency in Hz

option.WaitBar = false; % remove the wait bar of SIMUS
h = waitbar(0,'');
for k = 1:7
    waitbar(k/7,h,['SIMUS: RF series #' int2str(k) ' of 7'])
    RF{k} = simus(x,y,z,RC,txdel{k},param,option);
end
close(h)

% This is the 32th RF signal of the 1st series:

rf = RF{1}(:,32);
t = (0:numel(rf)-1)/param.fs*1e6; % time (ms)
plot(t,rf)
set(gca,'YColor','none','box','off')
xlabel('time (\mus)')
title('RF signal of the 32^{th} element (1^{st} series, tilt = -20{\circ})')
axis tight

%% Demodulate the RF signals with RF2IQ
% Before beamforming, the RF signals must be I/Q demodulated.

IQ = cell(7,1);  % this cell will contain the I/Q series

for k = 1:7
    IQ{k} = rf2iq(RF{k},param.fs,param.fc);
end

% This is the 32th I/Q signal of the 1st series:

iq = IQ{1}(:,32);
plot(t,real(iq),t,imag(iq))
set(gca,'YColor','none','box','off')
xlabel('time (\mus)')
title('I/Q signal of the 32^{th} element (1^{st} series, tilt = -20{\circ})')
legend({'in-phase','quadrature'})
axis tight

%% Beamform the I/Q signals with DAS

% To generate images of the left ventricle, beamform the I/Q signals onto a 256 $\times$ 128 polar grid.

% Generate the image (polar) grid using IMPOLGRID.

[xi,zi] = impolgrid([256 128],15e-2,deg2rad(80),param);

% Beamform the I/Q signals using a delay-and-sum with the function DAS.

bIQ = zeros(256,128,7);  % this array will contain the 7 I/Q images

h = waitbar(0,'');
for k = 1:7
    waitbar(k/7,h,['DAS: I/Q series #' int2str(k) ' of 7'])
    bIQ(:,:,k) = das(IQ{k},xi,zi,txdel{k},param);
end
close(h)

%% Time-gain compensate the beamformed I/Q with TGC

% Time-gain compensation tends to equalize the amplitudes along fast-time.

bIQ = tgc(bIQ);

%% Check the ultrasound images
% An ultrasound image is obtained by log-compressing the amplitude of the
% beamformed I/Q signals. Have a look at the images obtained when steering 
% at -20 degrees.

I = bmode(bIQ(:,:,1),50); % log-compressed image
pcolor(xi*1e2,zi*1e2,I)
shading interp, colormap gray
title('DW-based echo image with a tilt angle of -20{\circ}')

axis equal ij
set(gca,'XColor','none','box','off')
c = colorbar;
c.YTick = [0 255];
c.YTickLabel = {'-50 dB','0 dB'};
ylabel('[cm]')

% The individual images are of poor quality. The compound image obtained 
% with a series of 7 diverging waves steered at different angles is of 
% better quality:

cIQ = sum(bIQ,3); % this is the compound beamformed I/Q
I = bmode(cIQ,50); % log-compressed image
pcolor(xi*1e2,zi*1e2,I)
shading interp, colormap gray
title('Compound DW-based cardiac echo image')

axis equal ij
set(gca,'XColor','none','box','off')
c = colorbar;
c.YTick = [0 255];
c.YTickLabel = {'-50 dB','0 dB'};







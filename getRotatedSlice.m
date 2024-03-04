function rotatedSlice = getRotatedSlice(volume, angle, sliceIndex, axis)
    % Extracts a rotated slice from a 3D volume.
    %
    % Parameters:
    % - volume: 3D array, the CT volume.
    % - angle: double, the angle of rotation in degrees.
    % - sliceIndex: int, the index of the slice to extract after rotation.
    % - axis: char, the axis to rotate around ('x', 'y', or 'z').
    %
    % Returns:
    % - rotatedSlice: 2D array, the rotated slice.

    if nargin < 4
        axis = 'z'; % Default rotation axis is 'z' if not specified
    end

    % Validate input axis
    if ~ismember(axis, ['x', 'y', 'z'])
        error('Axis must be ''x'', ''y'', or ''z''.');
    end

    rotatedVolume = volume; % Initialize with original volume

    % Rotate the volume based on the specified axis
    switch axis
        case 'z'
            for i = 1:size(volume, 3)
                rotatedVolume(:, :, i) = imrotate(volume(:, :, i), angle, 'bilinear', 'crop');
            end
        case 'x'
            % For rotation around x, rotate each "vertical slice"
            for i = 1:size(volume, 1)
                rotatedVolume(i, :, :) = imrotate(squeeze(volume(i, :, :)), angle, 'bilinear', 'crop');
            end
        case 'y'
            % For rotation around y, rotate each "horizontal slice" from a side view
            for i = 1:size(volume, 2)
                rotatedVolume(:, i, :) = imrotate(squeeze(volume(:, i, :)), angle, 'bilinear', 'crop');
            end
    end

    % Extract the slice after rotation
    if strcmp(axis, 'z')
        rotatedSlice = rotatedVolume(:, :, sliceIndex);
    elseif strcmp(axis, 'x')
        rotatedSlice = squeeze(rotatedVolume(sliceIndex, :, :));
    elseif strcmp(axis, 'y')
        rotatedSlice = squeeze(rotatedVolume(:, sliceIndex, :));
    end
end

function dose = DeformDose(dose, vectorfield)
% DeformDose reads in a dose volume and 3D vector transformation matrix, 
% deforms the dose volume according to the matrix, and resamples the dose 
% volume back to the original (x,y,z) values using linear 3D interpolation.
% If a valid GPU device is available, the interpolation will be performed
% using CUDA.  Otherwise, DeformDose will revert to CPU interpolation.  In 
% place data manipulation is used where possible to improve execution 
% efficiency.
%
% The following variables are required for proper execution: 
%   dose: three-dimentional array containing the dose volume to be deformed
%   vectorfield: four-dimensional (x,y,z,3) array indicating the vector
%       transformations (dx,dy,dz) from the reference image to the merged
%       image.  The magnitude of each vector field is relative to the size
%       of one voxel in that dimension
%
% The following variables are returned upon succesful completion:
%   dose: three-dimensional array of deformed dose, resampled to
%       original coordinates
%
% Copyright (C) 2014 University of Wisconsin Board of Regents
%
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the  
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along 
% with this program. If not, see http://www.gnu.org/licenses/.

% If the size of the dose array does not equal the first three dimensions
% of the vector field
if size(dose,1) ~= size(vectorfield,1) || size(dose,2) ~= ...
        size(vectorfield,2) || size(dose,3) ~= size(vectorfield,3)
    
    % Report an error
    Event(['Error: the input dose array must be the same size as the ', ...
        'vector field'],'ERROR');
end

% Create unitless reference 3D meshgrids in x,y,z dimensions
[refX, refY, refZ] = meshgrid(1:size(vectorfield, 2), ...
    1:size(vectorfield, 1),1:size(vectorfield, 3));

% Log beginning of transformation and start timer
Event('Beginning deformed dose resampling using 3D linear interpolation');
tic;
    
% If the deformation field is null (in the case of deformation type NODIR,
% see DeformImages for more details)
if ~isequal(vectorfield, zeros(size(vectorfield)))
    % Attempt to use GPU; if it fails, catch the error and skip to CPU
    % method below
    try
        % Clear and initialize GPU memory
        gpuDevice(1);
        
        % Convert vector field to GPU array
        vectorfieldg = gpuArray(vectorfield);

        % Run interp3 to resample from modified meshgrids back to original
        % meshgrids, then gather GPU results back to CPU memory
        dose = gather(interp3(gpuArray(refX), gpuArray(refY), ...
            gpuArray(refZ), gpuArray(single(dose)), ...
            gpuArray(refX) + squeeze(vectorfieldg(:,:,:,1)), ...
            gpuArray(refY) + squeeze(vectorfieldg(:,:,:,2)), ...
            gpuArray(refZ) + squeeze(vectorfieldg(:,:,:,3)), 'linear', 0));
        
        % Log completion of interpolation
        Event(sprintf('GPU interpolation completed in %0.3f seconds', toc));
        
        % Clear GPU memory
        gpuDevice(1);
        
        % Clear temporary variable
        clear vectorfieldg;
        
    % Catch any errors when attempting to use GPU functions
    catch
        % Run interp3 to resample from modified meshgrids back to original
        % meshgrids
        dose = interp3(refX, refY, refZ, single(dose), ...
            refX + squeeze(vectorfield(:,:,:,1)), ...
            refY + squeeze(vectorfield(:,:,:,2)), ...
            refZ + squeeze(vectorfield(:,:,:,3)), 'linear', 0);
        
        % Log completion of CPU interpolation
        Event(sprintf('CPU interpolation completed in %0.3f seconds', toc), ...
            'WARN');
    end
end

% Clear temporary variables
clear refX refY refZ;

function vectorfield = DeformImages(referenceImage, mergedImage, method)
% DeformImages uses deformable image registration to deform a reference 
% image set to a daily/merged image set using various deformation methods
% and returns a four-dimensional vector transformation matrix.
%
% The following variables are required for proper execution: 
%   referenceImage: structure containing the image data, dimensions, width,
%       start coordinates, structure set UID, couch checksum and IVDT.  See
%       LoadReferenceImage for more information on contents of fields
%   mergedImage: structure containing the image data, dimensions, width,
%       start coordinates, IVDT, and dailyMask.  See MergeImages for more 
%       information
%   method: type of registration to perform.  See switch statement below
%       for options
%
% The following variables are returned upon succesful completion:
%   vectorfield: four-dimensional array of the same size as the reference
%       image indicating the vector transformations from the reference
%       image to the merged image.  The magnitude of each vector field is
%       relative to the size of one voxel in that dimension
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

% Debug flag.  If true, additional content will be output to log
debug = false;

% Log beginning of deformation and start timer
Event('Deforming reference image to merged daily image');
tic;

% Log which method was chosen
Event(['Method ', method, ' selected for deformation algorithm']);

% Convert reference image to equivalent merged image using linear
% interpolation of IVDTs
referenceImage.data = interp1(mergedImage.ivdt(:,2), mergedImage.ivdt(:,1), ...
    interp1(referenceImage.ivdt(:,1), referenceImage.ivdt(:,2), ...
    referenceImage.data, 'linear', 'extrap'), 'linear', 'extrap');

% Log density conversion
Event(['Reference image converted to merged-equivalent Hounsfield ', ...
    'Units using IVDT']);

% Based on the method passed to DeformImages, execute the code for that
% algorithm
switch method
    
%% CUDA-based Plastimatch B-Spline deformable registration based on MI
% This deformable registration technique requires plastimatch be installed
% on this workstation and compiled for GPU.  Data is passed to/from
% plastimatch using the ITK .mha file format and text command file.
% See http://iopscience.iop.org/0031-9155/55/21/001 for additional details
% on the algorithm.
case 'PM_BSPL_MI'
    %% Build reference MHA file
    % Generate a temprary filename for the reference image
    referenceFilename = [tempname, '.mha'];
    
    % Open a write file handle to the temporary reference image
    fid = fopen(referenceFilename, 'w', 'l');
    
    % Start writing the ITK header
    fprintf(fid, 'ObjectType=Image\n');
    fprintf(fid, 'NDims=3\n');
    
    % Specify the dimensions of the reference image
    fprintf(fid, 'DimSize=%i %i %i\n', referenceImage.dimensions);
    
    % Specify the data format (USHORT referring to unsigned 16-bit integer)
    fprintf(fid, 'ElementType=MET_USHORT\n');
    
    % Specify the byte order as little
    fprintf(fid, 'ElementByteOrderMSB=False\n');
    
    % Specify the reference voxel widths (in mm)
    fprintf(fid, 'ElementSize=%i %i %i\n', referenceImage.width*10);
    
    % Specify the reference voxel spacing to equal the widths (in mm)
    fprintf(fid, 'ElementSpacing=%i %i %i\n', referenceImage.width*10);
    
    % Specify the coordinate frame origin (in mm)
    fprintf(fid, 'Origin=%i %i %i\n', referenceImage.start*10);
    
    % Complete the .mha file header
    fprintf(fid, 'ElementDataFile=LOCAL\n');
    
    % Write the reference image data to the temporary file as uint16
    fwrite(fid, referenceImage.data, 'ushort', 0, 'l');
    
    % Close the file handle
    fclose(fid);
    
    % Clear the temporary variable
    clear fid;
    
    % Log where the reference file was saved
    Event(['Reference image written to ', referenceFilename]);
    
    %% Build merged MHA file
    % Generate a temporary file name for the merged image
    mergedFilename = [tempname, '.mha'];
    
    % Open a write file handle to the temporary merged image
    fid = fopen(mergedFilename, 'w', 'l');
    
    % Start writing the ITK header
    fprintf(fid, 'ObjectType=Image\n');
    fprintf(fid, 'NDims=3\n');
    
    % Specify the dimensions of the merged image
    fprintf(fid, 'DimSize=%i %i %i\n', mergedImage.dimensions);
    
    % Specify the data format (USHORT referring to unsigned 16-bit integer)
    fprintf(fid, 'ElementType=MET_USHORT\n');
    
    % Specify the byte order as little
    fprintf(fid, 'ElementByteOrderMSB=False\n');
    
    % Specify the merged voxel widths (in mm)
    fprintf(fid, 'ElementSize=%i %i %i\n', mergedImage.width*10);
    
    % Specify the merged voxel spacing to equal the widths (in mm)
    fprintf(fid, 'ElementSpacing=%i %i %i\n', mergedImage.width*10);
    
    % Specify the coordinate frame origin (in mm)
    fprintf(fid, 'Origin=%i %i %i\n', mergedImage.start*10);
    
    % Complete the .mha file header
    fprintf(fid, 'ElementDataFile=LOCAL\n');
    
    % Write the merged image data to the temporary file as uint16
    fwrite(fid, mergedImage.data, 'uint16', 0, 'l');
    
    % Close the file handle
    fclose(fid);
    
    % Clear the temporary variable
    clear fid;
    
    % Log where the merged file was saved
    Event(['Merged image written to ', mergedFilename]);
    
    %% Build merged MHA mask file
    % Generate a temporary file name for the merged image mask
    mergedMaskFilename = [tempname, '.mha'];
    
    % Open a write file handle to the temporary merged image mask
    fid = fopen(mergedMaskFilename, 'w', 'l');
    
    % Start writing the ITK header
    fprintf(fid, 'ObjectType=Image\n');
    fprintf(fid, 'NDims=3\n');
    
    % Specify the dimensions of the merged image
    fprintf(fid, 'DimSize=%i %i %i\n', mergedImage.dimensions);
    
    % Specify the data format (USHORT referring to unsigned 16-bit integer)
    fprintf(fid, 'ElementType=MET_USHORT\n');
    
    % Specify the byte order as little
    fprintf(fid, 'ElementByteOrderMSB=False\n');
    
    % Specify the merged voxel widths (in mm)
    fprintf(fid, 'ElementSize=%i %i %i\n', mergedImage.width*10);
    
    % Specify the merged voxel spacing to equal the widths (in mm)
    fprintf(fid, 'ElementSpacing=%i %i %i\n', mergedImage.width*10);
    
    % Specify the coordinate frame origin (in mm)
    fprintf(fid, 'Origin=%i %i %i\n', mergedImage.start*10);
    
    % Complete the .mha file header
    fprintf(fid, 'ElementDataFile=LOCAL\n');
    
    % Write the merged image mask data to the temporary file as uint16
    fwrite(fid, mergedImage.dailyMask, 'uint16', 0, 'l');
    
    % Close the file handle
    fclose(fid);
    
    % Clear the temporary variable
    clear fid;
    
    % Log where the merged file was saved
    Event(['Image mask written to ', mergedMaskFilename]);
    
    %% Build plastimatch command file
    % Generate a temporary file name for the command file
    commandFile = [tempname, '.txt'];
    
    % Open a write file handle to the temporary command file
    fid = fopen(commandFile, 'w');
    
    % Specify the inputs to the registration
    fprintf(fid, '[GLOBAL]\n');
    fprintf(fid, 'moving=%s\n', referenceFilename);
    fprintf(fid, 'fixed=%s\n', mergedFilename);
    fprintf(fid, 'moving_mask=%s\n', mergedMaskFilename);
    fprintf(fid, 'fixed_mask=%s\n', mergedMaskFilename);
    
    % If the debug flag is set
    if debug
        % Generate a temporary filename for the resulting image (for debug 
        % purposes) 
        outputimg = [tempname, '.mha']; %#ok<*UNRCH>

        % Specify the output image
        fprintf(fid, 'img_out=%s\n', outputimg);

        % Log deformed image location
        Event(['Deformed image will be saved to ', outputimg], 'DEBUG');

        % Clear temporary variables
        clear outputimg;
    end
    
    % Generate a temporary filename for the resulting coefficients 
    output = [tempname, '.txt'];
    
    % Specify the output file
    fprintf(fid, 'xform_out=%s\n', output);
    
    % Specify stage 1 deformable image registration parameters.  Refer to 
    % http://plastimatch.org/registration_command_file_reference.html for
    % more information on these parameters
    fprintf(fid, '[STAGE]\n');
    fprintf(fid, 'xform=align_center\n');
    
    % Specify stage 2 parameters
    fprintf(fid, '[STAGE]\n');
    fprintf(fid, 'impl=plastimatch\n');
    fprintf(fid, 'xform=bspline\n');
    fprintf(fid, 'optim=steepest\n');
    fprintf(fid, 'threading=cuda\n');
    fprintf(fid, 'metric=mi\n');
    fprintf(fid, 'max_its=100\n');
    fprintf(fid, 'grid_spac=30 30 30\n');
    fprintf(fid, 'regularization_lambda=0.000\n');
    fprintf(fid, 'res=4 4 2\n');
    
    % Specify stage 3 parameters
    fprintf(fid, '[STAGE]\n');
    fprintf(fid, 'impl=plastimatch\n');
    fprintf(fid, 'xform=bspline\n');
    fprintf(fid, 'optim=steepest\n');
    fprintf(fid, 'threading=cuda\n');
    fprintf(fid, 'metric=mi\n');
    fprintf(fid, 'max_its=100\n');
    fprintf(fid, 'grid_spac=20 20 20\n');
    fprintf(fid, 'regularization_lambda=0.000\n');
    fprintf(fid, 'res=2 2 1\n');
    
    % Specify stage 4 parameters
    fprintf(fid, '[STAGE]\n');
    fprintf(fid, 'impl=plastimatch\n');
    fprintf(fid, 'xform=bspline\n');
    fprintf(fid, 'optim=steepest\n');
    fprintf(fid, 'threading=cuda\n');
    fprintf(fid, 'metric=mi\n');
    fprintf(fid, 'max_its=50\n');
    fprintf(fid, 'grid_spac=10 10 10\n');
    fprintf(fid, 'regularization_lambda=0.000\n'); 
    fprintf(fid, 'res=1 1 1\n');
    
    % Close the file handle
    fclose(fid);
    
    % Clear temporary variables
    clear fid referenceFilename mergedFilename mergedMaskFilename;
    
    %% Run plastimatch
    % Log execution of system call
    Event(['Executing plastimatch register ', commandFile]);
    
    % Execute plastimatch using system call, saving the output and status
    [status, cmdout] = system(['plastimatch register ', commandFile]);
    
    % If the status == 0, the command completed successfully
    if status == 0
        % Log output
        Event(cmdout);
    else
        % Otherwise, plastimatch didn't complete succesfully, so log the 
        % resulting command output as an error
        Event(cmdout, 'ERROR');
    end
    
    % Clear temporary variables
    clear status cmdout commandFile;
    
    %% Convert transform into vector field
    % Generate a temporary filename for the vector field output file
    output2 = [tempname, '.mha'];
    
    % Log execution of system call
    Event(['Executing plastimatch xfconvert --input ', output, ...
        ' --output ', output2, ' --output-type vf']);
    
    % Execute plastimatch using system call, saving the output and status
    [status, cmdout] = system(['plastimatch xf-convert --input ', ...
        output, ' --output ', output2, ' --output-type vf']);
    
    % If the status == 0, the command completed successfully
    if status == 0
        % Log output
        Event(cmdout);
    else
        % Otherwise, plastimatch didn't complete succesfully, so log the 
        % resulting command output as an error
        Event(cmdout, 'ERROR');
    end
    
    % Clear temporary variables
    clear status cmdout output;
    
    %% Read in registration result  
    % Open file handle to temporary file
    fid = fopen(output2, 'r', 'l');
    
    % Skip over the ITK file header to place the file handle at the
    % beginning of the vector field data
    fseek(fid, -4 * prod(referenceImage.dimensions) * 3, 'eof');
    
    % Read in the vector field
    vectorfield = reshape(fread(fid, prod(referenceImage.dimensions) * 3, ...
        'single'), [3 referenceImage.dimensions(1) ...
        referenceImage.dimensions(2) referenceImage.dimensions(3)]);
    
    % Permute the vector field conents to put (x,y,z) first and the vector
    % magnitudes in the fourth dimension
    vectorfield = permute(vectorfield, [2 3 4 1]);
 
    % Log completion
    Event(sprintf(['Registration vectors read from %s with dimensions ', ...
        '(%i %i %i 3)'], output2, referenceImage.dimensions(1), ...
        referenceImage.dimensions(2), referenceImage.dimensions(3)));
    
    % Close the file handle
    fclose(fid);
    
    % Clear the temporary variables
    clear fid output2;
    
    % Normalize vectors to voxel units
    Event('Normalizing deformation vectors to unitless voxels');

    % Loop through each dimension
    for i = 1:3
        % Normalize the vector field dimension by the voxel size (in mm)
        vectorfield(:,:,:,i) = vectorfield(:,:,:,i) / ...
            referenceImage.width(i) / 10;
    end
    
%% Use Plastimatch DEMONS deformable registration using FSF
% This deformable registration technique requires plastimatch be installed
% on this workstation and computes via single-threaded CPU.  Data is passed 
% to/from plastimatch using the ITK .mha file format and text command file.
% See http://www.insight-journal.org/browse/publication/644 for additional
% detail on the algorithm.
case 'PM_DEMONS_FSF'
    %% Build reference MHA file
    % Generate a temprary filename for the reference image
    referenceFilename = [tempname, '.mha'];
    
    % Open a write file handle to the temporary reference image
    fid = fopen(referenceFilename, 'w', 'l');
    
    % Start writing the ITK header
    fprintf(fid, 'ObjectType=Image\n');
    fprintf(fid, 'NDims=3\n');
    
    % Specify the dimensions of the reference image
    fprintf(fid, 'DimSize=%i %i %i\n', referenceImage.dimensions);
    
    % Specify the data format (USHORT referring to unsigned 16-bit integer)
    fprintf(fid, 'ElementType=MET_USHORT\n');
    
    % Specify the byte order as little
    fprintf(fid, 'ElementByteOrderMSB=False\n');
    
    % Specify the reference voxel widths (in mm)
    fprintf(fid, 'ElementSize=%i %i %i\n', referenceImage.width*10);
    
    % Specify the reference voxel spacing to equal the widths (in mm)
    fprintf(fid, 'ElementSpacing=%i %i %i\n', referenceImage.width*10);
    
    % Specify the coordinate frame origin (in mm)
    fprintf(fid, 'Origin=%i %i %i\n', referenceImage.start*10);
    
    % Complete the .mha file header
    fprintf(fid, 'ElementDataFile=LOCAL\n');
    
    % Write the reference image data to the temporary file as uint16
    fwrite(fid, referenceImage.data, 'ushort', 0, 'l');
    
    % Close the file handle
    fclose(fid);
    
    % Clear the temporary variable
    clear fid;
    
    % Log where the reference file was saved
    Event(['Reference image written to ', referenceFilename]);
    
    %% Build merged MHA file
    % Generate a temporary file name for the merged image
    mergedFilename = [tempname, '.mha'];
    
    % Open a write file handle to the temporary merged image
    fid = fopen(mergedFilename, 'w', 'l');
    
    % Start writing the ITK header
    fprintf(fid, 'ObjectType=Image\n');
    fprintf(fid, 'NDims=3\n');
    
    % Specify the dimensions of the merged image
    fprintf(fid, 'DimSize=%i %i %i\n', mergedImage.dimensions);
    
    % Specify the data format (USHORT referring to unsigned 16-bit integer)
    fprintf(fid,'ElementType=MET_USHORT\n');
    
    % Specify the byte order as little
    fprintf(fid,'ElementByteOrderMSB=False\n');
    
    % Specify the merged voxel widths (in mm)
    fprintf(fid, 'ElementSize=%i %i %i\n', mergedImage.width*10);
    
    % Specify the merged voxel spacing to equal the widths (in mm)
    fprintf(fid, 'ElementSpacing=%i %i %i\n', mergedImage.width*10);
    
    % Specify the coordinate frame origin (in mm)
    fprintf(fid, 'Origin=%i %i %i\n', mergedImage.start*10);
    
    % Complete the .mha file header
    fprintf(fid, 'ElementDataFile=LOCAL\n');
    
    % Write the merged image data to the temporary file as uint16
    fwrite(fid, mergedImage.data, 'uint16', 0, 'l');
    
    % Close the file handle
    fclose(fid);
    
    % Clear the temporary variable
    clear fid;
    
    % Log where the merged file was saved
    Event(['Merged image written to ', mergedFilename]);
    
    %% Build merged MHA mask file
    % Generate a temporary file name for the merged image mask
    mergedMaskFilename = [tempname, '.mha'];
    
    % Open a write file handle to the temporary merged image mask
    fid = fopen(mergedMaskFilename, 'w', 'l');
    
    % Start writing the ITK header
    fprintf(fid, 'ObjectType=Image\n');
    fprintf(fid, 'NDims=3\n');
    
    % Specify the dimensions of the merged image mask
    fprintf(fid, 'DimSize=%i %i %i\n', mergedImage.dimensions);
    
    % Specify the data format (USHORT referring to unsigned 16-bit integer)
    fprintf(fid,'ElementType=MET_USHORT\n');
    
    % Specify the byte order as little
    fprintf(fid,'ElementByteOrderMSB=False\n');
    
    % Specify the merged voxel widths (in mm)
    fprintf(fid, 'ElementSize=%i %i %i\n', mergedImage.width*10);
    
    % Specify the merged voxel spacing to equal the widths (in mm)
    fprintf(fid, 'ElementSpacing=%i %i %i\n', mergedImage.width*10);
    
    % Specify the coordinate frame origin (in mm)
    fprintf(fid, 'Origin=%i %i %i\n', mergedImage.start*10);
    
    % Complete the .mha file header
    fprintf(fid, 'ElementDataFile=LOCAL\n');
    
    % Write the merged image mask data to the temporary file as uint16
    fwrite(fid, mergedImage.dailyMask, 'uint16', 0, 'l');
    
    % Close the file handle
    fclose(fid);
    
    % Clear the temporary variable
    clear fid;
    
    % Log where the merged mask file was saved
    Event(['Image mask written to ', mergedMaskFilename]);
    
    %% Build plastimatch command file
    % Generate a temporary file name for the command file
    commandFile = [tempname, '.txt'];
    
    % Open a write file handle to the temporary command file
    fid = fopen(commandFile, 'w');
    
    % Specify the inputs to the registration
    fprintf(fid, '[GLOBAL]\n');
    fprintf(fid, 'moving=%s\n', referenceFilename);
    fprintf(fid, 'fixed=%s\n', mergedFilename);
    fprintf(fid, 'moving_mask=%s\n', mergedMaskFilename);
    fprintf(fid, 'fixed_mask=%s\n', mergedMaskFilename);
    
    % Generate a temporary filename for the resulting coefficients 
    output = [tempname, '.mha'];
    
    % Specify the output file
    fprintf(fid, 'vf_out=%s\n', output);
    
    % Specify the image output type (single)
    fprintf(fid, 'img_out_type=short\n');
    
    % If the debug flag is set
    if debug
        % Generate a temporary filename for the resulting image (for debug 
        % purposes) 
        outputimg = [tempname, '.mha'];

        % Specify the output image
        fprintf(fid, 'img_out=%s\n', outputimg);

        % Log deformed image location
        Event(['Deformed image will be saved to ', outputimg], 'DEBUG');

        % Clear temporary variables
        clear outputimg;
    end
    
    % Specify stage 1 deformable image registration parameters.  Refer to 
    % http://plastimatch.org/registration_command_file_reference.html for
    % more information on these parameters
    fprintf(fid, '[STAGE]\n');
    fprintf(fid, 'xform=align_center\n');
    
    % Specify stage 2
    fprintf(fid, '[STAGE]\n');
    fprintf(fid, 'impl=plastimatch\n');
    fprintf(fid, 'xform=vf\n');
    fprintf(fid, 'optim=demons\n');
    fprintf(fid, 'optim_subtype=fsf\n');
    fprintf(fid, 'demons_gradient_type=symmetric\n');
    fprintf(fid, 'demons_smooth_deformation_field=1\n');
    fprintf(fid, 'demons_std_deformation_field=3\n');
    fprintf(fid, 'threading=openmp\n');
    fprintf(fid, 'max_its=200\n');
    fprintf(fid, 'res=4 4 2\n');
    
    % Specify stage 3
    fprintf(fid, '[STAGE]\n');
    fprintf(fid, 'impl=plastimatch\n');
    fprintf(fid, 'xform=vf\n');
    fprintf(fid, 'optim=demons\n');
    fprintf(fid, 'optim_subtype=fsf\n');
    fprintf(fid, 'demons_gradient_type=symmetric\n');
    fprintf(fid, 'demons_smooth_deformation_field=1\n');
    fprintf(fid, 'demons_std_deformation_field=2\n');
    fprintf(fid, 'threading=openmp\n');
    fprintf(fid, 'max_its=100\n');
    fprintf(fid, 'res=2 2 1\n');
    
    % Specify stage 4
    fprintf(fid, '[STAGE]\n');
    fprintf(fid, 'impl=plastimatch\n');
    fprintf(fid, 'xform=vf\n');
    fprintf(fid, 'optim=demons\n');
    fprintf(fid, 'optim_subtype=fsf\n');
    fprintf(fid, 'demons_gradient_type=symmetric\n');
    fprintf(fid, 'demons_smooth_deformation_field=0\n');
    fprintf(fid, 'demons_std_deformation_field=1\n');
    fprintf(fid, 'threading=openmp\n');
    fprintf(fid, 'max_its=100\n');
    fprintf(fid, 'res=1 1 1\n');
    
    %% Run plastimatch
    % Log execution of system call
    Event(['Executing plastimatch register ', commandFile]);
    
    % Execute plastimatch using system call, saving the output and status
    [status, cmdout] = system(['plastimatch register ', commandFile]);
    
    % If the status == 0, the command completed successfully
    if status == 0
        % Log output
        Event(cmdout);
    else
        % Otherwise, plastimatch didn't complete succesfully, so log the 
        % resulting command output as an error
        Event(cmdout, 'ERROR');
    end
    
    % Clear temporary variables
    clear status cmdout commandFile;
    
    %% Read in registration result  
    % Open file handle to temporary file
    fid = fopen(output, 'r', 'l');
    
    % Skip over the ITK file header to place the file handle at the
    % beginning of the vector field data
    fseek(fid, -4 * prod(referenceImage.dimensions) * 3, 'eof');
    
    % Read in the vector field
    vectorfield = reshape(fread(fid, prod(referenceImage.dimensions) * 3, ...
        'single'), [3 referenceImage.dimensions(1) ...
        referenceImage.dimensions(2) referenceImage.dimensions(3)]);
    
    % Permute the vector field conents to put (x,y,z) first and the vector
    % magnitudes in the fourth dimension
    vectorfield = permute(vectorfield, [2 3 4 1]);
 
    % Log completion
    Event(sprintf(['Registration vectors read from %s with dimensions', ...
        ' (%i %i %i 3)'], ...
        output, referenceImage.dimensions));
    
    % Close the file handle
    fclose(fid);
    
    % Clear the temporary variables
    clear fid output;
    
    % Normalize vectors to voxel units
    Event('Normalizing deformation vectors to unitless voxels');

    % Loop through each dimension
    for i = 1:3
        % Normalize the vector field dimension by the voxel size (in mm)
        vectorfield(:,:,:,i) = vectorfield(:,:,:,i) / ...
            referenceImage.width(i) / 10;
    end
    
%% No deformation
% This algorithm is used in demonstration mode when the system does not
% have any deformation libraries installed
case 'NODIR'
    % Return empty vector field
    vectorfield = zeros(referenceImage.dimensions(1), ...
        referenceImage.dimensions(2), referenceImage.dimensions(3), 3);
    
% Otherwise, an unsupported method was passed
otherwise
    % Log error
    Event(['Unsupported method ', method, ' passed to DeformImages'], ...
        'ERROR');
end

% Log completion of deformable image registration
Event(sprintf('Deformable image registration completed in %0.3f seconds', ...
    toc));

% Clear temporary variables
clear referenceImage mergedImage method debug;
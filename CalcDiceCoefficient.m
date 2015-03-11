function varargout = CalcDiceCoefficient(fixed, moving)
% CalcDiceCoefficient computes the Dice coefficient and average Hausdorff
% boundary distance for a set of structures given a fixed and moving mask.
% The statistics are computed using plastimatch.  For more information, see
% http://plastimatch.org/plastimatch.html#plastimatch-dice.
%
% The location of the moving data set is computed by adding the difference
% between the registration and rigid translations to the start coordinates.
% Accounting for rotations is currently not supported, and will generate a
% warning with 'N/A' returned in the return variables.
%
% The following variables are required for proper execution: 
%   fixed: a structure containing the following fields for the reference  
%       (or fixed) data: name, dimensions, start (in cm), width (in cm), 
%       and structures, where structures{i}.mask is a 3D logical mask of 
%       the contour's voxels.  For more information, see Loadfixed and 
%       LoadReferenceStructures.
%   moving: a structure containing the following fields for the
%       merged/daily (or moving) data: dimensions, start (in cm), width (in 
%       cm), registration, rigid, and structures, where structures{i}.mask 
%       is a 3D logical mask of the contour's voxels.  For more 
%       information, see MergeImage and DeformStructures.
%
% The following variables are returned upon succesful completion:
%   varargout{1}: if nargout == 1, an array of Dice coefficients for each
%       structure in the fixed and moving cell arrays
%   varargout{2}: if nargout == 2, an array of Hausdorff distances for each
%       structure in the fixed and moving cell arrays, in mm
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

% Initialize error flag.  If set to true, 'N/A' will be returned
na = false;

% Check for plastimatch
[~,cmdout] = system('which plastimatch');

% If plastimatch was not found
if strcmp(cmdout,'')   
    % Set the error flag to true
    na = true;
    
    % Log error
    Event(['Contour statistics could not be computed because ', ...
        'plastimatch is not available'], 'WARN');
end

% Clear temporary variable
clear cmdout;

% Check if any rotation/rigid adjustments are present
if ~isfield(moving,'registration') || ~isfield(moving,'rigid') || ...
        moving.registration(1) ~= moving.registration(1) || ...
        moving.registration(2) ~= moving.registration(2) || ...
        moving.registration(3) ~= moving.registration(3)
    % Set the error flag to true
    na = true;
    
    % Log error
    Event(['Contour statistics could not be computed because non-zero ', ...
        'registration adjustments are not currently supported'], 'WARN');
end

% Initialize return variables
dice = cell(1, size(fixed.structures, 2));
hausdorff = cell(1, size(fixed.structures, 2));

% If the error flag is still false, continue with computations
if ~na
    Event('Computing contour dice coefficients');
    tic;
    
    % Compute registration difference (in cm)
    adj = zeros(1,3);
    adj(1) = moving.registration(4) - moving.rigid(4);
    adj(2) = moving.registration(6) - moving.rigid(6);
    adj(3) = moving.registration(5) - moving.rigid(5);
    
    % Loop through each structure
    for i = 1:size(fixed.structures, 2)
        %% Build fixed MHA file
        % Generate a temprary filename for the fixed structure mask
        referenceFilename = [tempname, '.mha'];

        % Open a write file handle to the temporary fixed image
        fid = fopen(referenceFilename, 'w', 'l');

        % Start writing the ITK header
        fprintf(fid, 'ObjectType=Image\n');
        fprintf(fid, 'NDims=3\n');

        % Specify the dimensions of the fixed image
        fprintf(fid, 'DimSize=%i %i %i\n', fixed.dimensions);

        % Specify the data format (USHORT referring to unsigned 16-bit integer)
        fprintf(fid, 'ElementType=MET_USHORT\n');

        % Specify the byte order as little
        fprintf(fid, 'ElementByteOrderMSB=False\n');

        % Specify the fixed voxel widths (in mm)
        fprintf(fid, 'ElementSize=%i %i %i\n', fixed.width * 10);

        % Specify the fixed voxel spacing to equal the widths (in mm)
        fprintf(fid, 'ElementSpacing=%i %i %i\n', fixed.width * 10);
    
        % Specify the coordinate frame origin (in mm)
        fprintf(fid, 'Origin=%i %i %i\n', fixed.start * 10);

        % Complete the .mha file header
        fprintf(fid, 'ElementDataFile=LOCAL\n');

        % Write the structure mask to the temporary file as uint16
        fwrite(fid, fixed.structures{i}.mask, 'ushort', 0, 'l');

        % Close the file handle
        fclose(fid);

        % Clear the temporary variable
        clear fid;

        %% Build moving MHA file
        % Generate a temprary filename for the moving structure mask
        mergedFilename = [tempname, '.mha'];

        % Open a write file handle to the temporary moving image
        fid = fopen(mergedFilename, 'w', 'l');

        % Start writing the ITK header
        fprintf(fid, 'ObjectType=Image\n');
        fprintf(fid, 'NDims=3\n');

        % Specify the dimensions of the moving image
        fprintf(fid, 'DimSize=%i %i %i\n', moving.dimensions);

        % Specify the data format (USHORT referring to unsigned 16-bit integer)
        fprintf(fid, 'ElementType=MET_USHORT\n');

        % Specify the byte order as little
        fprintf(fid, 'ElementByteOrderMSB=False\n');

        % Specify the moving voxel widths (in mm)
        fprintf(fid, 'ElementSize=%i %i %i\n', moving.width * 10);

        % Specify the moving voxel spacing to equal the widths (in mm)
        fprintf(fid, 'ElementSpacing=%i %i %i\n', moving.width * 10);
    
        % Specify the coordinate frame origin (in mm), applying
        % difference between registration and rigid adjustment
        fprintf(fid, 'Origin=%i %i %i\n', (moving.start + adj) * 10);

        % Complete the .mha file header
        fprintf(fid, 'ElementDataFile=LOCAL\n');

        % Write the structure mask to the temporary file as uint16
        fwrite(fid, moving.structures{i}.mask, 'ushort', 0, 'l');

        % Close the file handle
        fclose(fid);

        % Clear the temporary variable
        clear fid;
    
        % If only one return variable is required, only compute DICE
        if nargout == 1
            % Execute plastimatch using system call, saving the output
            [status, cmdout] = system(['plastimatch dice --dice ', ...
                referenceFilename, ' ', mergedFilename]);
            
            % Parse Dice coefficient
            [~, tokens] = regexp(cmdout, 'DICE:[ ]+([0-9\.]+)', 'match', 'tokens');
            
            % If the status == 0, the command completed successfully
            if status == 0 && size(tokens,1) > 0
                % Set return variables
                dice{i} = sprintf('%0.3f', str2double(tokens{1}));
            else
                % Return 'N/A'
                dice{i} = 'N/A';
            end
            
            % Clear temporary variable
            clear tokens;
            
            % Log result
            Event(['Structure ', fixed.structures{i}.name, ...
                ' Dice coefficient computed as ', dice{i}]);
        
        % Otherwise, compute both DICE and Hausdorff distance
        else
            % Execute plastimatch using system call, saving the output
            [status, cmdout] = system(['plastimatch dice --dice ', ...
                referenceFilename, ' ', mergedFilename]);
            
            % Parse Dice coefficient
            [~, tokens] = regexp(cmdout, 'DICE:[ ]+([0-9\.]+)', 'match', 'tokens');
            
            % If the status == 0, the command completed successfully
            if status == 0 && size(tokens,1) > 0
                % Set return variables
                dice{i} = sprintf('%0.3f', str2double(tokens{1}));
            else
                % Return 'N/A'
                dice{i} = 'N/A';
            end
            
            % Clear temporary variable
            clear tokens;
            
            % Execute plastimatch using system call, saving the output
            [status, cmdout] = system(['plastimatch dice --hausdorff ', ...
                referenceFilename, ' ', mergedFilename]);
            
            % Parse Dice coefficient
            [~, tokens] = regexp(cmdout, 'Average Hausdorff distance \(boundary\) = ([0-9\.]+)', 'match', 'tokens');
            
            % If the status == 0, the command completed successfully
            if status == 0 && size(tokens,1) > 0
                % Set return variables
                hausdorff{i} = sprintf('%0.3f', str2double(tokens{1}));
            else
                % Return 'N/A'
                hausdorff{i} = 'N/A';
            end
            
            % Clear temporary variable
            clear tokens;
            
            % Log result
            Event(['Structure ', fixed.structures{i}.name, ...
                ' Dice coefficient and Hausdorff average distance computed as ', ...
                dice{i}, ' and ', hausdorff{i}]);
        end
        
        % Clear temporary variables
        clear status cmdout;
    end
    
    % Clear temporary variables
    clear i adj;
    
% Otherwise, the error flag was set to true
else 
    % Return 'N/A' for each value
    for i = 1:size(fixed.structures, 2)
        dice{i} = 'N/A';
        hausdorff{i} = 'N/A';
    end
    
    % Clear temproary variables
    clear i;
end
   
% Log completion
Event(sprintf('Dice coefficient computation completed succesfully in %0.3f seconds', toc));

% Set varargout
if nargout >= 1; varargout{1} = dice; end
if nargout >= 2; varargout{2} = hausdorff; end

% Clear temporary variables
clear na dice hausdorff;
    
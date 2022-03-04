% ***********************************************************************************
%             I S O M A P   V I S U A L I Z A T I O N   E X A M P L E
% ***********************************************************************************
%
%  This is the script to execute an IsoMap dimension reducing visualization of two
%  different sets of images from the train_set folder. The two different sets are
%  face and nonface images and the images are read as a set of 19x19 greyscale
%  values and then each 19x19 matrix representing a single image is reshaped to be 
%  one 361 element which are stacked to create the X matrix. The IsoMap algorithm 
%  is then performed on this X matrix. Most of this script is just loading the 
%  images and constructing the X matrix to perform IsoMap on. The end result is a 
%  3-D visulaization where each face and nonface image is represented by a point
%  with blue representing the face images and red indicating the nonface images
%
%  This was originally written by Giovanni Nuno and can be found at:
% 
%                      https://github.com/gionuno/isomap
%
%  Here we have commented the code so it is clear how this code can be adapted to 
%  a distance preserving IsoMap application to patient specific HFpEF and HFrEF 
%  optimized parameter values.
%
%  Code commented by:      Brian Carlson
%                          Physiolgical Systems Dynamics Lab
%                          Department of Molecular and Integrative Physiology
%                          Univrsity of Michigan
%
%  Initially created on:   3 March 2022
%  Modified on:            4 March 2022
% 
% ***********************************************************************************

    % Grabs the path to the folders in 'train_set' and then initializes matrices
    DIR = dir('train_set');
    X = [];
    C = [];

    % Since the . and .. paths are also loaded into DIR we just want the last two
    %  which in this example are 'face' and 'nonface'
    for i = 3:length(DIR)

        % Builds the path to the 'face' and 'nonface' folders
        dirname = strcat('train_set/',DIR(i).name);
        % Moves into the 'face' or 'nonface' folder and loads all the
        %  names of the files in that folder
        AUX_DIR = dir(dirname);
        % Again ignoring the . and .. paths this steps through
        %  all the images in the current folder
        for j = 3:length(AUX_DIR)

            % Makes the path to image j-3 and reads it as a 19x19
            %  array of greyscale values. The greyscale value is divided
            %  by 255 to put it on the scale of 0 to 1
            auxfilename = strcat(dirname,'/',AUX_DIR(j).name);
            aux = double(imread(auxfilename))/255.0;
            % Reshape the 19x19 matrix as a row vector and then add it to the
            %  matrix X which has rows representing the images and the columns
            %  representing the normalized greyscale value of the regions in the image
            X = [X ; reshape(aux,[1,numel(aux)])];
            % Assigning the color blue to all face and red to all nonface images
            if (i == 3)
                C = [C ; [0 0 256]];
            else
                C = [C ; [256 0 0]];
            end

        end

    end

    % Calling IsoMap to operate on X by first making a Kd-tree, searching for the
    %  nearest neighbor for each point in X with the tree
    %  to preserve distance and dimension reduced to d
    K = 32;                    % Number of nearest neighbors to use in knnsearch
    d = 3;                     % Number of dimensions to reduce to
    [Y,idxNN,D] = isomap(X,K,d);

    % Visualize the IsoMap transformation
    scatter3(Y(:,1),Y(:,2),Y(:,3),d^2,C);
    
    
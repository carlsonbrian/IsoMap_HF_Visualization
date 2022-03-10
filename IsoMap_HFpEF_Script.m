% ***********************************************************************************
%            I S O M A P   H F p E F / H F r E F   V I S U A L I Z A T I O N
% ***********************************************************************************
%
%  This is the script uses the IsoMap dimension reducing algorithm as implemented 
%  by Giovanni Nuno to visualize our HFpEF/HFrEF patient-specific optimized 
%  parameter space in 2-D. We will then be able to see if this distance preserving
%  representation is similar to what we get with PCA telling us that the underlying
%  data structure is adequately approximated with a linear approach.


%  The two different sets are
%  face and nonface images and the images are read as a set of 19x19 greyscale
%  values and then each 19x19 matrix representing a single image is reshaped to be 
%  one 361 element which are stacked to create the X matrix. The IsoMap algorithm 
%  is then performed on this X matrix. Most of this script is just loading the 
%  images and constructing the X matrix to perform IsoMap on. The end result is a 
%  3-D visulaization where each face and nonface image is represented by a point
%  with blue representing the face images and red indicating the nonface images
%
%  This code was adapted from the original code by Giovanni Nuno which can be 
%  found at:
% 
%                      https://github.com/gionuno/isomap
%
%  Code adapted by:        Brian Carlson and Edith Jones
%                          Physiolgical Systems Dynamics Lab
%                          Department of Molecular and Integrative Physiology
%                          Univrsity of Michigan
%
%  Initially created on:   9 March 2022
%  Last modified on:       9 March 2022
% 
% ***********************************************************************************

    %% Load optimized parameter values

    % Load optimized parameters from text file
    %  69 patients represented in rows, with the first 9 columns being
    %  the optimized parameters, the 10 column being the HF type and the
    %  last column being the patient number
    load HFpEFvsHFrEF_Optp.txt 
    
    % Split out the parameters from HF type and patient number
    A_Optp = HFpEFvsHFrEF_Optp(:,1:9);              % Optim parameters
    HFType = HFpEFvsHFrEF_Optp(:,10);               % Heart failure type
    PatNum = HFpEFvsHFrEF_Optp(:,11);               % Patient number
    
    % Finding out the total number of HFrEF and HFpEF patients
    Num_Pats = size(A_Optp,1);                      % Number of patients
    Num_PatsHFrEF = 0;                              % Initialize HFrEF counter
    C = zeros(Num_Pats,3);                          % Init plotting color vector
    for i = 1:Num_Pats                              % Count up all patients
        if (HFType(i) == 0)                         %  with HF type 0 (HFrEF)
            Num_PatsHFrEF = Num_PatsHFrEF + 1;
            C(i,1) = 256;                           % HFrEF are red
        end
    end
    Num_PatsHFpEF = 0;                              % Initialize HFpEF counter
    for i = 1:Num_Pats                              % Count up all patients
        if (HFType(i) == 1)                         %  with HF type 1 (HFrEF)
            Num_PatsHFpEF = Num_PatsHFpEF + 1;
            C(i,3) = 256;                           % HFrEF are blue
        end
    end
    
    % Grabbing the number of optimized parameters and assigning
    %  the parameter name to each column for plotting purposes
    Num_Optp = size(A_Optp,2);                      % Number of optimized parameters
    Optp_Names = {'E_{LV}', '\lambda_{LV}', ...     % Parameter names
        'E_{RV}','\lambda_{RV}','E_{PA}', ...
        'E_{PV}','R_{pul}','E_{SA}','R_{sys}'};

    
%% Normalizing data

    % This normailzation flag can be changed to reflect how we are normalizing
    %  this dataset
    %  Norm_Flag = 0 - no normalization 
    %            = 1 - normalize by mean
    %            = 2 - normalize by standard deviation
    Norm_Flag = 2;    
    
    % Taking the mean of each optimized parameter and centering the matrix
    AMean_Optp = mean(A_Optp,1);                % Measure mean across all patients
    ACent_Optp = zeros(Num_Pats, Num_Optp);     % Preallocate
    for i = 1:Num_Pats
        for j = 1:Num_Optp
            ACent_Optp(i,j) = A_Optp(i,j) - ... % Centering the data matrix
                AMean_Optp(j);    
        end
    end
    
    % Now normalizing (or not) depending on normalization flag
    AStd_Optp = std(A_Optp);                    % Measure std dev across all pats
    ANorm_Optp = zeros(Num_Pats, Num_Optp);     % Preallocate
    if (Norm_Flag == 0)                         % No normalization
        ANorm_Optp = ACent_Optp;
    elseif (Norm_Flag == 1)                     % Normalized by mean
        for i = 1:Num_Pats
            for j = 1:Num_Optp
                ANorm_Optp(i,j) = ...
                    ACent_Optp(i,j) / AMean_Optp(j);
            end
        end
    else                                        % Normalized by standard deviation
        for i = 1:Num_Pats
            for j = 1:Num_Optp
                ANorm_Optp(i,j) = ...
                    ACent_Optp(i,j) / AStd_Optp(j);
            end
        end
    end
    
    
    %% IsoMap code
    
    % Setting the number of nearest neighbors and the dimension to reduce to.
    %  A good rule of thumb is to set the nearest neighbors to sqrt(NumPats) so
    %  initially we have set K equal to sqrt(69) ~ 8
    K = 8;                      % Number of nearest neighbors to use in knnsearch
    d = 2;                      % Number of dimensions to reduce to
    
    % Creating the k-dimensional tree from the ANorm_Optp matrix containing 
    %  the patients as rows and the normalized and centered optimized parameters
    %  as columns (if normalized, if not the optimized parametes are just centered)
    Tree = createns(ANorm_Optp,'NSMethod','kdtree');
    
    % Now using this kd tree model we can search over it to find the nearest
    %  neighbor of each patient in ANorm_Optp. This decreases the time to find 
    %  the nearest neighbor since we are throwing out most of the patients that 
    %  are not close to our target patient on the kd tree. We are adding an 
    %  additional nearest neighbor to the list (K+1) since the first column 
    %  will contain the target point as nearest to itself
    idxNN = knnsearch(Tree,ANorm_Optp,'K',K+1);
    % This just cuts off the first column leaving us with an list of
    %  indices of the K nearest neighbors to each patient (row) in ANorm_Optp
    idxNN = idxNN(:,2:end);
    
    % Preallocates a matrix W that is an Num_PatsxNum_Pats matrix full of 
    %  infinities. When the W matrix is filled up with euclidean distances 
    %  from the K nearest neighbors below it is assumed the distance to all
    %  other patients is infinite to begin with.
    W = inf*ones(Num_Pats,Num_Pats);
    % This steps through all N images and calculates the L2 norm
    %  between image i and its K nearest neighbors
    for i = 1:Num_Pats
        for k = 1:K
            W(i,idxNN(i,k)) = norm(ANorm_Optp(i,:) - ...
                ANorm_Optp(idxNN(i,k),:));
        end
    end
    
    % This Floyd-Warshall algorithm recursively steps through the distance
    %  matrix W and replaces any direct distance between image i and image j with 
    %  a distance through another image k if the distance from i to k + k to j
    %  is smaller. Here it is important that not too big or too small of a number
    %  of nearest neighbors is selected. Too big and we may find a shorter distance
    %  that is not supported by the inherent structure of the dataset. Too small
    %  and we introduce noise to our solution. Note that this FW algortihm is 
    %  performed on min(W,W') which just says that the starting matrix should be 
    %  symmetric with the distance from i to j equal to the distance from j to i 
    %  (e.g. no directionality enforced)
    D = floyd_warshall(min(W,W'));
    
    % We use J to center the element-wise squared distance matrix and then
    %  apply an eigenvalue decomposition on the result to find the first d
    %  eigenvectors with the largest eigenvalues to represent the IsoMap in
    %  euclidean d space.
    J = eye(Num_Pats)-ones(Num_Pats,Num_Pats)/Num_Pats;
    [Y,~] = eig(-0.5*J*(D.^2)*J);
    
    % Visualize the IsoMap transformation
    if (d == 2)
        scatter(Y(:,1),Y(:,2),d^3,C);
    elseif (d == 3)
        scatter3(Y(:,1),Y(:,2),Y(:,3),d^2,C);
    end

 
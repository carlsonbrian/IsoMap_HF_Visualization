%% ***********************************************************************************
%                         I S O M A P   F U N C T I O N
% ***********************************************************************************
%  In this function is an IsoMap dimension reducing transformation of the data 
%  utilizing two existing Matlab functions createns and knnsearch. In createns a
%  kd tree model is built which partitions the dataset into regions that are 
%  rouhgly equally populated based on the original dimensionality of the dataset. 
%  This model structure can then be used to search for the K nearest neighbors 
%  using the knnsearch function. With the K nearest neighbors of each point we can 
%  then construct a directed graph weighted by the shortest geodesic distance 
%  between each point using the Floyd-Warshall algorithm. With this final geodesic
%  distance matrix an eigenvalue decomposition is then used to map back into a 
%  euclidean space of dimension d for visualization.
%
%  This was originally written by Giovanni Nuno and can be found at:
%                   https://github.com/gionuno/isomap
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
%% ***********************************************************************************

function [Y,idxNN,D] = isomap(X,K,d)

    % Creating the k-dimensional tree from the X matrix containing the 
    %  images as rows and the greyscale values of discretized regions of
    %  the images as columns
    Tree = createns(X,'NSMethod','kdtree');

    % Now using this kd tree model we can search over it to find the nearest
    %  neighbor of each image in X. This decreases the time to find the nearest
    %  neighbor since we are throwing out most of the images that are not close
    %  to our target image on the kd tree. We are adding an additional nearest
    %  neighbor to the list since the first column will contain the target point
    %  as nearest to itself
    idxNN = knnsearch(Tree,X,'K',K+1);
    % This just cuts off the first column leaving us with an list of
    %  indices of the K nearest neighbors to each image (row) in X
    idxNN = idxNN(:,2:end);

    % Finds the number of images in X and then preallocates a
    %  matrix W that is an NxN matrix full of infinities. When the
    %  W matrix is filled up with euclidean distances from the K
    %  nearest neighbors below it is assumed the distance to all
    %  other images is infinite to begin with.
    N = size(X,1);
    W = inf*ones(N,N);

    % This steps through all N images and calculates the L2 norm
    %  between image i and its K nearest neighbors
    for i = 1:N
        disp(i);
        for k = 1:K
            W(i,idxNN(i,k)) = norm(X(i,:)-X(idxNN(i,k),:));
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
    %  euclidean d space. (Not sure why the first eigenvector is skipped here)
    J = eye(N)-ones(N,N)/N;
    [Y,~] = eig(-0.5*J*(D.^2)*J);
    Y = Y(:,2:d+1);
    
end
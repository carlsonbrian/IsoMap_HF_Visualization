%% ***********************************************************************************
%                    F L O Y D - W A R S H A L L   F U N C T I O N
% ***********************************************************************************
%  In this function a matrix containing the euclidean distance of each point to
%  its K nearest neighbors is the input and using these distances a shortest path
%  from each point to all other points is constructed. This new matrix containing
%  the geodesic distances is then output for further transformation into the
%  reduced dimension for visualization.
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

    function D = floyd_warshall(D_0)
    
        % Start with the nearest neighbor matrix where all points outside the 
        %  nearest neighbor are given a path length of infinity
        D = D_0;
        N = size(D,1);
        
        % Simple loop which fills up the distances to the points not a K
        %  nearest neighbor to point k with the geodesic distance through 
        %  other points in the dataset. This does this by stepping through
        %  each point one by one and updating D_0 to D_1 by looking at all
        %  paths that go through point 1 and seeing if they are shorter than
        %  the distances already in D_0, etc. The repmat construct below does
        %  this in a non-intuitive way
        for k = 1:N
            D = min(D,repmat(D(:,k),[1,N])+repmat(D(k,:),[N,1]));
        end
        
    end
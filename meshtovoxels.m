function mask = meshtovoxels(varargin)
    % meshtovoxels Voxelize a closed mesh
    %
    %   Syntax
    %   ------
    %   mask = meshvolume('f', F, 'v', V, 'x', X, 'y', Y, 'z', Z)
    %
    %   Description
    %   -----------
    %   mask = meshvolume('f', F, 'v', V, 'x', X, 'y', Y, 'z', Z) generates
    %   a boolean mask from the mesh represented by faces F and vertices V
    %   using the grid vectors X, Y and Z. The faces input F is an mx3 or
    %   3xm array that specifies the vertices to connect in the vertices
    %   array V. X, Y and Z represent the grid vectors to use to create the
    %   meshgrid for voxel creation.
    %
    %   Adapted from: Fast parallel surface and solid voxelization on GPUs,
    %   Schwarz and Seidel, ACM SIGGRAPH ASIA 2010 PAPERS
    %   DOI: 10.1145/1882261.1866201
    %
    %   © 2016, Peter Beemiller (pbeemiller@gmail.com)
    %
    %   See also SurfacesReader | isosurface | patch
    
    %% Parse the inputs.
    parser = inputParser;
    parser.KeepUnmatched = true;
    
    parser.addParameter('faces', [], ...
        @(f)ismatrix(f) && any(size(f) == 3))
    parser.addParameter('vertices', [], ...
        @(v)ismatrix(v) && any(size(v) == 3))
    parser.addParameter('xgridvector', [], @(X)isvector(X))
    parser.addParameter('ygridvector', [], @(Y)isvector(Y))
    parser.addParameter('zgridvector', [], @(Z)isvector(Z))
    
    parser.parse(varargin{:})
    
    %% Get the faces and vertices inputs.
    faces = parser.Results.faces;
    vertices = parser.Results.vertices;
    
    %% Convert the arrays to mx3 if necessary.
    if size(faces, 1) == 3
        faces = transpose(faces);
    end % if
    
    if size(vertices, 1) == 3
        vertices = transpose(vertices);
    end % if
    
    %% Create the triangles array and calculate the bounding boxes.
    triangles(:, :, 1) = vertices(faces(:, 1), :);
    triangles(:, :, 2) = vertices(faces(:, 2), :);
    triangles(:, :, 3) = vertices(faces(:, 3), :);
    
    triMin = min(triangles, [], 3);
    triMax = max(triangles, [], 3);
    
    %% Calculate the mesh bounding box.
    meshMin.x = min(vertices(:, 1));
    meshMin.y = min(vertices(:, 2));
    meshMin.z = min(vertices(:, 3));

    meshMax.x = max(vertices(:, 1));
    meshMax.y = max(vertices(:, 2));
    meshMax.z = max(vertices(:, 3));
    
    %% Create the grids.
    xVector = parser.Results.xgridvector;
    yVector = parser.Results.ygridvector;
    zVector = parser.Results.zgridvector;
    
    dp = [...
        xVector(2) - xVector(1), ...
        yVector(2) - yVector(1), ...
        zVector(2) - zVector(1)];
    
    [xGrid, yGrid, zGrid] = meshgrid(xVector, yVector, zVector);
    
    %% Calculate the mesh bounding box voxel mask.
    box = ...
        xGrid >= meshMin.x - dp(1) & xGrid <= meshMax.x + dp(1) & ...
        yGrid >= meshMin.y - dp(2) & yGrid <= meshMax.y + dp(2) & ...
        zGrid >= meshMin.z - dp(3) & zGrid <= meshMax.z + dp(3);
    
    % Only test voxels in the bounding box.
    idxBox = find(box);
    
    %% Allocate the mask.
    mask = false(size(box));
    
    %% Fill voxels that enclose a vertex.
    for v = idxBox'
        % Calculate the voxel bounding box.
        p = [xGrid(v), yGrid(v), zGrid(v)];
        pM = p + dp;
        
        if any(...
            (triMax(:, 1) - p(1)) .* (triMin(:, 1) - pM(1)) < 0.0 & ...
            (triMax(:, 2) - p(2)) .* (triMin(:, 2) - pM(2)) < 0.0 & ...
            (triMax(:, 3) - p(3)) .* (triMin(:, 3) - pM(3)) < 0.0)
            
            mask(v) = true;
        end % if
    end % for v
    
    %% Remove the voxels that have been filled from the search.
    idxBox = setdiff(idxBox, find(mask));
    
    %% Fill voxels that are traversed by an edge.
    for t = 1:size(triangles, 1)
        % Get the triangle vertices.
        v1 = triangles(t, :, 1);
        v2 = triangles(t, :, 2);
        v3 = triangles(t, :, 3);
        
        % Calculate the triangles edge vectors in voxel coordinates.
        e1 = v2 - v1;
        e2 = v3 - v2;
        e3 = v1 - v3;
        
        % Calculate the triangle plane normal.
        n = cross(e1, e2);
        sN = sign(n);
        sN(sN == 0) = 1;

        %% Find voxels that overlap the triangle plane.
        % Get the voxel coordinates.
        pV = [xGrid(idxBox), yGrid(idxBox), zGrid(idxBox)];
        
        % Get the test points.
        c = dp;
        c(n <= 0) = 0;
        
        % Calculate the offsets.
        d1 = dot(n, c - v1);
        d2 = dot(n, (dp - c) - v1);
        
        dnp = dot(repmat(n, [size(pV, 1), 1]), pV, 2);
        
        isPlane = (dnp + d1).*(dnp + d2) <= 0;
        
        idxPlane = idxBox(isPlane);
        
        if isempty(idxPlane)
            continue
        end % if
        
        %% Find voxels that overlap the 2D projections.
        pV = [xGrid(idxPlane), yGrid(idxPlane), zGrid(idxPlane)];
        
        xyOverlap = ...
            (triMax(t, 1) - pV(:, 1)) .* (triMin(t, 1) - pV(:, 1) - dp(1)) < 0 & ...
            (triMax(t, 2) - pV(:, 2)) .* (triMin(t, 2) - pV(:, 2) - dp(2)) < 0;
        
        yzOverlap = ...
            (triMax(t, 3) - pV(:, 2)) .* (triMin(t, 2) - pV(:, 2) - dp(3)) < 0 & ...
            (triMax(t, 3) - pV(:, 3)) .* (triMin(t, 3) - pV(:, 3) - dp(3)) < 0;
        
        xzOverlap = ...
            (triMax(t, 1) - pV(:, 1)) .* (triMin(t, 1) - pV(:, 1) - dp(1)) < 0 & ...
            (triMax(t, 3) - pV(:, 3)) .* (triMin(t, 3) - pV(:, 3) - dp(3)) < 0;
        
        idxOverlap = idxPlane(xyOverlap & yzOverlap & xzOverlap);
        
        if isempty(idxOverlap)
            continue
        end % if
        
        %% Setup edge voxel intersection test.
        % XY plane:
        xyEN = sN(3)*[...
            -e1(2), e1(1);
            -e2(2), e2(1);
            -e3(2), e3(1)];
        
        xyD = -dot(xyEN, [v1(1:2); v2(1:2); v3(1:2)], 2) + ...
            max(cat(2, [0; 0; 0], dp(1)*xyEN(:, 1)), [], 2) + ...
            max(cat(2, [0; 0; 0], dp(2)*xyEN(:, 2)), [], 2);
        
        % YZ plane:
        yzEN = sN(1)*[...
            -e1(3), e1(2);
            -e2(3), e2(2);
            -e3(3), e3(2)];

        yzD = -dot(yzEN, [v1(2:3); v2(2:3); v3(2:3)], 2) + ...
            max(cat(2, [0; 0; 0], dp(2)*yzEN(:, 1)), [], 2) + ...
            max(cat(2, [0; 0; 0], dp(3)*yzEN(:, 2)), [], 2);

        %xz plane:
        xzEN = -sN(2)*[...
            -e1(3), e1(1);
            -e2(3), e2(1);
            -e3(3), e3(1)];

        xzD = -dot(xzEN, [v1([1, 3]); v2([1, 3]); v3([1, 3])], 2) + ...
            max(cat(2, [0; 0; 0], dp(1)*xzEN(:, 1)), [], 2) + ...
            max(cat(2, [0; 0; 0], dp(3)*xzEN(:, 2)), [], 2);

        %% Check for overlap of the 2D projections.
        for v = idxOverlap'
            p = [xGrid(v), yGrid(v), zGrid(v)];
            
            if edgeoverlap('xy') && edgeoverlap('yz') && edgeoverlap('xz')
                mask(v) = true;
                idxBox(idxBox == v) = [];
                break
            end % if
        end % for v
    end % for t

    %% Nested function for edge overlap tests
    function tf = edgeoverlap(plane)
        % overlap2d Test triangle edge 2D projections for voxel overlap
        %
        %
        
        %% Calculate the edge normals in the desired plane.
        switch plane
            
            case 'xy'
                tf = all(dot(xyEN, [p(1:2); p(1:2); p(1:2)], 2) + xyD > 0);

            case 'yz'
                tf = all(dot(yzEN, [p(2:3); p(2:3); p(2:3)], 2) + yzD > 0);

            case 'xz'
                tf = all(dot(xzEN, [p([1, 3]); p([1, 3]); p([1, 3])], 2) + xzD > 0);

        end % switch
    end % edgeoverlap    
end % meshtovoxels

function generate_dam_break_h5_mesh(factor)

nx = 10; ny = 5;

nx = nx*factor; ny = ny*factor;

% number of vertices in x and y dir
nvx = nx + 1;
nvy = ny + 1;

% 1D array for the location of vertices
xv_1d = [0:nx]/factor;
yv_1d = [0:ny]/factor;

% 2D arrays for vertex coordinates
[xv_2d, yv_2d] = meshgrid(xv_1d,yv_1d);
xv_2d = xv_2d';
yv_2d = yv_2d';

% create a base mesh consisting of quad cells in which all cells are active
base_vert_id       = reshape([1:nvx*nvy],nvx,nvy);
base_ncells        = nx*ny;
max_faces_per_cell = 4;
base_cells         = zeros(base_ncells,max_faces_per_cell);

count = 0;
for jj = 1:ny
    for ii = 1:nx
        count = count + 1;
        base_cells(count,:) = [...
            base_vert_id(ii  ,jj  ) ...
            base_vert_id(ii+1,jj  ) ...
            base_vert_id(ii+1,jj+1) ...
            base_vert_id(ii  ,jj+1)
            ];
    end
end

base_xc = mean(xv_2d(base_cells),2);
base_yc = mean(yv_2d(base_cells),2);

% create a mask for active/inactive cells
mask = ones(nx,ny);
% loc = find(base_xc > 4 & base_xc < 6 & base_yc > 4 & base_yc < 5); mask(loc) = 0;
% loc = find(base_xc > 4 & base_xc < 6 & base_yc > 0 & base_yc < 2); mask(loc) = 0;

if (nx == 5 && ny == 4)
    mask(2,[1 3 4]) = 0;
    mask(3,[1 3 4]) = 0;
end

% find active cells
loc_valid_cells = find(mask == 1);
ncells = length(loc_valid_cells);
fprintf('ncells    = %d\n',ncells);

% map the cells from base to new mesh
cell_id_base2new = mask*0;
cell_id_base2new(loc_valid_cells) = [1:ncells];

if (0)
    figure;
    plot(xv_2d,yv_2d,'ok','markerfacecolor','k','markersize',10)
    hold all
    for icell = 1:base_ncells
        idx = [base_cells(icell,:) base_cells(icell,1)];
        plot(xv_2d(idx), yv_2d(idx), '-b', 'linewidth', 2);
    end

    for icell = 1:ncells
        idx = [base_cells(loc_valid_cells(icell),:) base_cells(loc_valid_cells(icell),1)];
        fill(xv_2d(idx), yv_2d(idx), 'r');
    end
end

% let's identify the vertices that are active
vertex_mask = base_vert_id*0;

% set 1 in the vertex mask for vertices that are part of active cells
for ivert = 1:max_faces_per_cell
    id = base_cells(loc_valid_cells,ivert);
    vertex_mask(id) = 1;
end

% find the total number of vertices that are active
nvertices = sum(vertex_mask(:));
fprintf('nvertices = %d\n',nvertices);

% map the vertex IDs from the base mesh to new mesh
loc_valid_vertices = find(vertex_mask == 1);
vertex_id_base2new = vertex_mask*0;
vertex_id_base2new(loc_valid_vertices) = [1:nvertices];

% for each active cells in the new mesh, get the new vertex ID
cells_new = vertex_id_base2new(base_cells(loc_valid_cells,:));

% get the coordinates of active vertices
xv_new = xv_2d(loc_valid_vertices);
yv_new = yv_2d(loc_valid_vertices);

if (0)
    figure
    plot(xv_new,yv_new,'ok','markerfacecolor','k','markersize',10)
    hold all
    for icell = 1:ncells
        idx = [cells_new(icell,:) cells_new(icell,1)];
        plot(xv_new(idx), yv_new(idx), '-b', 'linewidth', 2);
    end
end

%                                                 DMPlex number
%
%  v6 --- f8 --- v7 --- f11--- v8       10 --- 21 --- 11 --- 24 --- 12
%  |             |             |        |             |             |
%  |             |             |        |             |             |
%  f9     c2     f7     c3     f10      22     2      20      3     23
%  |             |             |        |             |             |
%  |             |             |        |             |             |
%  v3 --- f2 --- v4 --- f6 --- v5  ===> 7  --- 15 --- 8  --- 19 --- 9
%  |             |             |        |             |             |
%  |             |             |        |             |             |
%  f3     c0     f1     c1     f5       16     0      14     1      18
%  |             |             |        |             |             |
%  |             |             |        |             |             |
%  v0 --- f0 --- v1--- f4 --- v2        4  --- 13 --- 5 --- 17 ---  6
%
%
% 9  vertices
% 12 faces
% 4  cells


% Let's number the faces in the order shown on above the left schematic i.e. in `f*`
%
%  c0:  f0 f1  f2   f3
%  c1:  f4 f5  f6  -f1
%  c2: -f2 f7  f8   f9
%  c3: -f6 f10 f11 -f7
%
%  Note: only the first or the fourth faces could be negative.
%

new_face_id_onBaseMesh     = base_cells*0;
face_count      = 0;
base_cell_count = 0;

% on the new mesh: for each face, keep track of vertex IDs that form the
% face
face_to_vertex = zeros(base_ncells*max_faces_per_cell, 2);

disp('finding faces');
for jj = 1:ny
    for ii = 1:nx
        base_cell_count = base_cell_count + 1;

        % only consider base cell that is active
        if (mask(ii,jj))

            %
            % 1st face of the cell
            %

            if (jj == 1) % Is the cell southmost?                
                % add a new face
                face_count = face_count + 1;             
                new_face_id_onBaseMesh(base_cell_count, 1) = face_count;

                % let's find the new vertex IDs forming the face

                %  - first, get the new cell ID
                id = cell_id_base2new(ii,jj);

                %  - next, get the new vertex IDs
                face_to_vertex(face_count,:) = cells_new(id,[1 2]);

                %cells_new_to_face_new(id, 1) = face_count;

            else

                % Check if the cell south of the current cell is active
                if (mask(ii,jj-1) == 1)

                    % Since the cell south of the current cell is active,
                    % use the 3rd face ID of the cell south of the current
                    % cell.
                    new_face_id_onBaseMesh(base_cell_count, 1) = -new_face_id_onBaseMesh(base_cell_count-nx, 3);

                else
                    % The cell south of the current cell was inactive, so
                    % add a new face.

                    % add a new face
                    face_count = face_count + 1;
                    new_face_id_onBaseMesh(base_cell_count, 1) = face_count;

                    % find new vertex IDs for the cell
                    id = cell_id_base2new(ii,jj);
                    face_to_vertex(face_count,:) = cells_new(id,[1 2]);
                end
            end
            
            %
            % Add 2nd face of the cell
            %
            face_count = face_count + 1;
            new_face_id_onBaseMesh(base_cell_count, 2) = face_count;
            id = cell_id_base2new(ii,jj);
            face_to_vertex(face_count,:) = cells_new(id,2:3);

            %
            % Add 3rd face of the cell
            %
            face_count = face_count + 1;
            new_face_id_onBaseMesh(base_cell_count, 3) = face_count;
            id = cell_id_base2new(ii,jj);
            face_to_vertex(face_count,:) = cells_new(id,3:4);

        
            % 4th face
            if (ii == 1) % Is the cell westmost cell?

                % Add a new face
                face_count = face_count + 1;
                new_face_id_onBaseMesh(base_cell_count, 4) = face_count;

                id = cell_id_base2new(ii,jj);
                face_to_vertex(face_count,:) = cells_new(id,[4 1]);

            else
                % Check if the cell west of the current cell is active
                if (mask(ii-1,jj) == 1)

                    % Since the cell west of the current cell is active,
                    % use the 2nd face of the cell west of the current cell
                    new_face_id_onBaseMesh(base_cell_count, 4) = -new_face_id_onBaseMesh(base_cell_count-1, 2);
                else
                    % Add a new face
                    face_count = face_count + 1;
                    new_face_id_onBaseMesh(base_cell_count, 4) = face_count;
                    id = cell_id_base2new(ii,jj);
                    face_to_vertex(face_count,:) = cells_new(id,[4 1]);
                end
            end
        end
    end
end

nfaces = max(max(new_face_id_onBaseMesh));
fprintf('nfaces    = %d\n',nfaces);

offset_vertices = ncells;
offset_faces    = offset_vertices + nvertices;

% extract new face IDs for active cells
new_face_id = new_face_id_onBaseMesh(loc_valid_cells,:);

% put the face IDs in 1D format
new_face_id_1d = reshape(abs(new_face_id_onBaseMesh(loc_valid_cells,:))', ncells*max_faces_per_cell,1);

% put the face-to-vertex in 1D format
face_to_vetex_1d = reshape(face_to_vertex(1:nfaces,:)',nfaces*2,1);

% Now put the data in the DMPlex format

% 1. Create `cells`
cells = [
    new_face_id_1d-1+offset_faces;      % subtract 1 to convert from 1-based to 0-based indexing; add offset for face
    face_to_vetex_1d-1+offset_vertices; % subtract 1 to convert from 1-based to 0-based indexing; add offset for vertices
];

% 2. Create `cone`
cones = zeros(ncells + nfaces + nvertices,1);
cones(1                     :ncells            ) = max_faces_per_cell;  % cells   : number of faces
cones(ncells + 1            :ncells + nvertices) = 0;                   % vertices: zero cones for the vertices
cones(ncells + nvertices + 1:end               ) = 2;                   % faces   : two vertices forming the face

% 3. `order` 
order = [0:ncells + nfaces + nvertices-1]';

% 4. `orientation`
orientation = cells*0;
loc = find(new_face_id' < 0);
orientation(loc) = -1;       % set -1 for the faces

% 5. vertex coordinates
vertices = [xv_new yv_new]';

out_fname = sprintf('dam_break_%dx%d.h5',nx,ny);
disp(out_fname);

write_h5mesh_file(out_fname, vertices, cells, cones, order, orientation);

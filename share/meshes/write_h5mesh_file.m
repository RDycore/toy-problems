function write_h5mesh_file(out_fname, vertices, cells, cones, order, orientation, ncells)

system(sprintf('rm -rf %s',out_fname));

h5create(out_fname,'/topology/cells',[1 length(cells)],'Datatype','int32');
h5create(out_fname,'/topology/cones',[1 length(cones)],'Datatype','int32');
h5create(out_fname,'/topology/order',[1 length(order)],'Datatype','int32');
h5create(out_fname,'/topology/orientation',[1 length(orientation)],'Datatype','int32');

h5create(out_fname,'/labels/Cell Sets/1/indices',[1 ncells],'Datatype','int32');

h5writeatt(out_fname,'/topology/cells','cell_dim',[size(vertices,1)])

h5create(out_fname,'/geometry/vertices',size(vertices));

h5write(out_fname,'/topology/cells',int32(cells'))
h5write(out_fname,'/topology/cones',int32(cones'))
h5write(out_fname,'/topology/order',int32(order'))
h5write(out_fname,'/topology/orientation',int32(orientation'))
h5write(out_fname,'/geometry/vertices',vertices);

h5write(out_fname,'/labels/Cell Sets/1/indices',int32([0:ncells-1]));


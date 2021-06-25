function h5writeSpatialRef(filename,location,rasterref)
% h5writeSpatialRef(filename,location,rasterref)
%
% filename - HDF 5 file
% location - e.g. Group appropriate to this spatial reference
% rasterref - MATLAB spatial reference object

% just the relevant fields of the input object
RasterSize = rasterref.RasterSize;
assert(contains(class(rasterref),'rasterref'),...
    'rasterref is class %s, not recognized',class(rasterref))
RefMatrix = RasterRef2RefMat(rasterref);

h5writeatt(filename,location,'RasterClass',class(rasterref));
h5writeatt(filename,location,'RasterSize',RasterSize);
h5writeatt(filename,location,'RefMatrix',RefMatrix);

end
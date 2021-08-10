function model = loadmodel(fname)
model.x     = h5read(fname,'/Terrace/x');
model.z     = h5read(fname,'/Terrace/z');
model.U     = h5read(fname,'/U0');
model.betaz	= h5read(fname,'/beta_z');
model.slope = h5read(fname,'/slope');
end
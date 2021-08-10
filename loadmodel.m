function model = loadmodel(fname)
model.x     = h5read(fname,'/Terrace/x');
model.z     = h5read(fname,'/Terrace/z');
model.U     = h5read(fname,'/U0');
model.betaz	= h5read(fname,'/beta_z');
model.slope = h5read(fname,'/slope');
<<<<<<< HEAD
end
=======
end
>>>>>>> 1a08b648b484c24b7bb8751b7b6220f259ac72c1

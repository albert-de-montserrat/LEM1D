% function model = loadmodel(fname)
% A           = h5read(fname,1);
% model.x     = A(:,1);
% model.z     = A(:,2);
% model.U     = A(1,4);
% model.betaz	= A(1,5);
% model.hwb   = A(1,6);
% model.slope = A(1,7);
% end
function model = loadmodel(fname)
model.x     = h5read(fname,'/Terrace/x');
model.z     = h5read(fname,'/Terrace/z');
model.zr    = h5read(fname,'/Rivers/z');
% model.age   = h5read(fname,'/Age/age');
% model.re    = h5read(fname,'/Time/reoccupation_time');
% model.U     = h5read(fname,'/Parameters/U0');
% model.betaz	= h5read(fname,'/Parameters');
% model.slope = h5read(fname,'/Parameters/slope');
end
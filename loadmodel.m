function model = loadmodel(fname)
A           = csvread(fname,1);
model.x     = A(:,1);
model.z     = A(:,2);
model.U     = A(1,4);
model.betaz	= A(1,5);
model.hwb   = A(1,6);
model.slope = A(1,7);
end
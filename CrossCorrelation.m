% -- Load reference model
[xref,zref]     = loadreference('OmuraProfile.mat');
[xref,isort]    = sort(xref);
zref            = zref(isort);

% -- Load Julia CSV model
model = loadmodel('potato.h5');
hsea  = 101.33;

% -- Resample reference model
dx    = model.x(2)-model.x(1);
xnew  = xref(1):dx:xref(end);
zref  = interp1(xref,zref,xnew); xref=xnew;

% -- Crop Julias model
ikeep   = model.z <= 250 & model.z >= -200;
model.x = model.x(ikeep); model.x = model.x - model.x(1);
model.z = model.z(ikeep); model.z = model.z - hsea;

% -- Cross correlation
[r,lags] = xcorr(zref,model.z);

% -- Plots 
figure(1);clf;hold on
plot(xref,zref,'r')
plot(model.x,model.z,'bx')

figure(2);
plot(lags,r,'r')

function [xref,zref] = loadreference(fname)
A       = load(fname);
xref    = A.x;
zref    = A.z;
end

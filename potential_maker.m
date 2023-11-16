clear; close all; clc;

%a = 2.84Å. see const.m for more stuff
a1=[const.a,0];
a2=[0,const.a];
a3=[0,0,const.a];
[b1,b2,b3] = Reciprocal([a1,0],[a2,0],a3);
%Number of grid points, number of Z points, and number of lattices
%contained in the overall superlattice (or rather the square root of that)
Ncell = 32; Nz = 100; Nslat = 2;
zMax = 6; zMin = -2;%units Å

V = zeros(Ncell,Ncell,Nz);
X = zeros(Ncell,Ncell);
Y = zeros(Ncell,Ncell);
Z = linspace(zMin,zMax,Nz);
for i = 1:Ncell
    for j = 1:Ncell
        X(i,j) = (a1(1)*i+a2(1)*j)./Ncell;
        Y(i,j) = (a1(2)*i+a2(2)*j)./Ncell;
    end
end

for z = 1:Nz
    V(:,:,z) = Vfunc(X,Y,Z(z));
end

%We strictly ought to be careful with boundary conditions cos MS doesn't
%actually check them lol
%===

Vsuper = zeros(Nslat*Ncell,Nslat*Ncell,Nz);
for z = 1:Nz
    for nx = 1:Ncell:Nslat*Ncell
        for ny = 1:Ncell:Nslat*Ncell
            Vsuper(nx:nx+Ncell-1,ny:ny+Ncell-1,z) = V(:,:,z);
        end
    end
end
writematrix(Vsuper,"V.csv")

potStructArray.V = Vsuper;

Multiscat.PreparePotentialFiles(potStructArray);

Multiscat.prepareFourierLabels(Vsuper);

potStructArray.a1=a1; potStructArray.a2=a2;
potStructArray.zmin=Z(1);
potStructArray.zmax=Z(end);
potStructArray.zPoints=length(Z);


confStruct=Multiscat.createConfigStruct(potStructArray);
Multiscat.prepareConfigFile(confStruct);

%We also prepare a .csv which contains an equipotential plot.
%===
equipotValue = 0;%Units meV ig
equipotentialMat = zeros(Ncell*Nslat,Ncell*Nslat);
M = max(Vsuper,[],"all");
for i = 1:Ncell*Nslat
    for j = 1:Ncell*Nslat
        foundVal = false;
        for k = 1:Nz
            if(Vsuper(i,j,k) < equipotValue && ~foundVal)
                disp("i, j, k = " + i + ", " + j + ", " + k +...
                    ", Value = " + Vsuper(i,j,k) +...
                    ", Z(k) = " + Z(k));
                equipotentialMat(i,j) = Z(k);
                foundVal = true;
            end
        end
        if(~foundVal)
            %error("unable to find equipotential!")
        end
    end
end
%disp(equipotentialMat)
writematrix(equipotentialMat,'Equipotential.csv','Delimiter', ',')
%===

function [b1,b2,b3] = Reciprocal(a1,a2,a3)
    factor = 2*pi/dot(cross(a1,a2),a3);
    b1 = factor*cross(a2,a3);
    b2 = factor*cross(a3,a1);
    b3 = factor*cross(a1,a2);
end

function [VmatrixElement] = Vfunc(x,y,z)
    D = 7.63;
    alpha = 1.1;
    z0 = 1.0;
    beta = 0.1;
    function [V0] = V0func(z)
        V0 = D * exp(2*alpha*(z0-z))...
            - 2*D*exp(alpha*(z0-z));
    end
    function [V1] = V1func(z)
        V1 = -2*beta*D*exp(2*alpha*(z0-z));
    end
    function [Q] = Qfunc(x,y)
        Q = cos(2*pi*x/const.a) + cos(2*pi*y/const.a);
    end
    VmatrixElement = V0func(z) + V1func(z)...
        * Qfunc(x,y);
end
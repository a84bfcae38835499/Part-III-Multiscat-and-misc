clear; close all; clc;

%a = 2.84Å. see const.m for more stuff
%a1=[const.a,0];
%a2=[0,const.a];
a1=[const.b,0];
a2=[const.b/2,const.b * sqrt(3)/2];
a3=[0,0,const.b];
[b1,b2,b3] = Reciprocal([a1,0],[a2,0],a3);
%Number of grid points, number of Z points, and number of lattices
%contained in the overall superlattice (or rather the square root of that)
Ncell = 64; Nz = 100; Nsuper = 2;
zMax = 6; zMin = -2;%units Å

disp("M * [0,1] = ")
disp(const.sheerMat*[0;1])

V = zeros(Ncell,Ncell,Nz);
X = zeros(Ncell,Ncell);
Y = zeros(Ncell,Ncell);
Xsuper = zeros(Ncell*Nsuper,Ncell*Nsuper);
Ysuper = zeros(Ncell*Nsuper,Ncell*Nsuper);
Z = linspace(zMin,zMax,Nz);

for i = 1:Ncell
    for j = 1:Ncell
        X(i,j) = (a1(1)*i+a2(1)*j)./Ncell;
        Y(i,j) = (a1(2)*i+a2(2)*j)./Ncell;
    end
end

for i = 1:Ncell*Nsuper
    for j = 1:Ncell*Nsuper
        Xsuper(i,j) = (a1(1)*i+a2(1)*j)./Ncell;
        Ysuper(i,j) = (a1(2)*i+a2(2)*j)./Ncell;
    end
end
  for k = 1:Nz
        V(:,:,k) = Vfunc(X,Y,Z(k));
  end
%We strictly ought to be careful with boundary conditions cos MS doesn't
%actually check them lol
%===
%% Now we duplicate the lattice to get a superlattice

Vsuper = zeros(Nsuper*Ncell,Nsuper*Ncell,Nz);
for z = 1:Nz
    for nx = 1:Ncell:Nsuper*Ncell
        for ny = 1:Ncell:Nsuper*Ncell
            Vsuper(nx:nx+Ncell-1,ny:ny+Ncell-1,z) = V(:,:,z);
        end
    end
end
%writematrix(Vsuper,"V.csv")

%Vsuper = readmatrix("V_boyao.csv");
%Vsuper = reshape(Vsuper,[Ncell,Ncell,Nz]);

%===
%% Now add imperfections to the lattice
if false
  for k = 1:size(V,3) %Should be the z direction
    dropoff = Dropoff(Z(k));
    for i = 1:Ncell*Nsuper
      for j = 1:Ncell*Nsuper
        x = Xsuper(i,j);
        y = Ysuper(i,j);
        val = Gaussian2D(x,y, ...
          [const.a*1.3,const.a*1.3],const.a/2,3*const.D*dropoff);
        Vsuper(i,j,k) = Vsuper(i,j,k)+val;
        disp("x, y, z = " + x + ", " + y + ", " + Z(k) +...
            ", Value = " + val);
      end
    end
  end
end
%===
%% We also prepare a .csv which contains an equipotential plot.
equipotValue = 0;%Units meV ig
eqCharArr = [num2str(equipotValue,'%+g') , ' meV'];
equipotentialMat = zeros(Ncell*Nsuper,Ncell*Nsuper);
M = max(Vsuper,[],"all");
for i = 1:Ncell*Nsuper
    for j = 1:Ncell*Nsuper
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

%===
%% All this stuff is just to write a comment line lmao
writematrix(equipotentialMat,'Equipotential.csv','Delimiter', ',')
S = fileread('Equipotential.csv');
S = ['#Energy surface at ', eqCharArr, newline, S];
FID = fopen('Equipotential.csv', 'w');
if FID == -1, error('Cannot open file %s', FileName); end
fwrite(FID, S, 'char');
fclose(FID);

%===
%% Plot the potential
% Plot of a slice of the potential in the nth row, that is for constant x
row = floor(Ncell/2);
figure
contourf(Z,  linspace(0, sqrt(3)*const.b*Nsuper, Ncell*Nsuper), ...%!!!
    reshape(Vsuper(row,:,:), [Ncell*Nsuper,Nz]), linspace(-20,100,24))

    fontsize(gcf,scale=1)
xlabel('z/Å')
ylabel('y/Å') %is this x or y? I think y but idrk
colorbar
xlim([-1,4])
title('Potential in z, used in simulation')
hbar = colorbar;
ylabel(hbar,'Energy / meV');

%% Plot the potential
% Linearly interpolated equipotential plot
    fontsize(gcf,scale=1)
equipotential_plot('V', Vsuper, 'z', Z, 'X', Xsuper, 'Y', Ysuper)

%% Plot the potential
zSample = 1.0;
zRow = floor((zSample - zMin)/(zMax-zMin) * Nz)
figure
contourf(Xsuper,Ysuper,Vsuper(:,:,zRow),10)
daspect([1 1 1])
xlabel('x/Å')
ylabel('y/Å')
title('Potentials at z = ' + string(zSample) + ' Å');
colormap(parula(10))
hbar = colorbar('southoutside');
xlabel(hbar,'Energy / meV');

    fontsize(gcf,scale=1)
%===
%% We supply the lattice to the mulitscat script so it can do its thing

    potStructArray.V = Vsuper;

Multiscat.PreparePotentialFiles(potStructArray);

Multiscat.prepareFourierLabels(Vsuper);

potStructArray.a1=a1; potStructArray.a2=a2;
potStructArray.zmin=Z(1);
potStructArray.zmax=Z(end);
potStructArray.zPoints=length(Z);

confStruct=Multiscat.createConfigStruct(potStructArray);
Multiscat.prepareConfigFile(confStruct);

%===
%% Function definitions

function [b1,b2,b3] = Reciprocal(a1,a2,a3)
    factor = 2*pi/dot(cross(a1,a2),a3);
    b1 = factor*cross(a2,a3);
    b2 = factor*cross(a3,a1);
    b3 = factor*cross(a1,a2);
end

function [VmatrixElement] = Vfunc(X,Y,Z)
    function [V0] = V0func(z)
        V0 = const.D * exp(2*const.alpha*(const.z0-z))...
            -2*const.D*exp(const.alpha*(const.z0-z));
    end
    function [V1] = V1func(z)
        V1 = -2*const.beta*const.D*exp(2*const.alpha*(const.z0-z));
    end
    function [V2] = V2func(z)
        V2 = -2*const.beta*const.D*exp(2*const.alpha*(const.z0-z));
    end
    function [Q] = Qfunc(x,y)
        Q = cos(2*pi*x/const.b) + cos(2*pi*y/const.b);
    end
  function [Q] = QhexfuncSingle(x,y)
        %disp("[][][][][]")
        %disp(v1)
        %disp(const.sheerMat)
        %disp("[][][][][]")
        %nu = y * 2/sqrt(3);
        %mu = x - (y/(sqrt(3)));
        x_n = x / (3*const.b);
        y_n = y / (sqrt(3)*const.b);
        Q = 0;
        
        mu_n1 = x_n*2;
        nu_n1 = y_n - x_n;
        Q = Q + cos(2*pi*nu_n1) + cos(2*pi*mu_n1);

        mu_n2 = x_n*2;
        nu_n2 = -y_n - x_n;
        Q = Q + cos(2*pi*nu_n2) + cos(2*pi*mu_n2);

        nu_n3 = y_n - x_n;
        mu_n3 = -y_n - x_n;
        Q = Q + cos(2*pi*nu_n3) + cos(2*pi*mu_n3);
        Q = Q/3;
        %Q = cos(2*pi*nu/const.a)^5 + cos(2*pi*mu/const.a)^5;
    end
  function [Q] = Qhexfunc(X,Y)
        X_n = X ./ (const.b);
        Y_n = Y ./ (const.b/sqrt(3));
        Q = zeros(length(X),length(Y));
        
        mu = X_n.*2;
        nu = Y_n - X_n;
        Q = Q + cos(2*pi*nu) + cos(2*pi*mu);

        %mu = X_n.*2;
        nu = -Y_n - X_n;
        Q = Q + cos(2*pi*nu) + cos(2*pi*mu);

        nu = Y_n - X_n;
        mu = -Y_n - X_n;
        Q = Q + cos(2*pi*nu) + cos(2*pi*mu);
        
        Q = Q./3;
        %Q = cos(2*pi*nu/const.a)^5 + cos(2*pi*mu/const.a)^5;
    end
    VmatrixElement = V0func(Z) ...
        + V1func(Z) * Qhexfunc(X,Y)...
        + Qhexfunc(X-const.b/3,Y) * V2func(Z);
end

function [DV] = Dropoff(z)
  %Use this to attenuate the gaussian in z
    DV = -exp(2*const.alpha*(const.z0-z));
end

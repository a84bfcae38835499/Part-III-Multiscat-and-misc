clear; close all; clc;
rng default;

%a = 2.84Å. see const.m for more stuff
%a1=[const.a,0];
%a2=[0,const.a]; 
a1=[const.b,0]; %These lattice parameters correspond to the projected hexagon tiling of MoS2, which is neither the bond length (because that's in 3D)
                %nor the unit cell vectors, which are three times these
                %lengths
a2=[const.b/2,const.b * sqrt(3)/2];
a3=[0,0,const.b];
%A1 = a1;
%A2 = a2;
A1 = a1 * 3;
A2 = a2 * 3;
[b1,b2,b3] = Reciprocal([a1,0],[a2,0],a3);
[B1,B2,B3] = Reciprocal([A1,0],[A2,0],a3);
%% data for python hex plotter WIP
writematrix([],'latticeVects.info_for_vivian_python_nice_plotting_hexagon_script',FileType='text')
A1str = [char(num2str(A1))];
A2str = [char(num2str(A2))];
B1str = [char(num2str(B1(1:2)))];
B2str = [char(num2str(B2(1:2)))];
S = fileread('latticeVects.info_for_vivian_python_nice_plotting_hexagon_script');
realStr = ['Real space vectors:',newline,'A1 = ',A1str, newline, 'A2 = ',A2str];
recpStr = ['Reciprocal vectors:',newline,'B1 = ',B1str, newline, 'B2 = ', B2str];

S = [realStr,newline,recpStr,S];
FID = fopen('latticeVects.info_for_vivian_python_nice_plotting_hexagon_script', 'w');
if FID == -1, error('Cannot open file %s', FileName); end
fwrite(FID, S, 'char');
fclose(FID);
%%


%Number of grid points, number of Z points, and number of lattices
%contained in the overall superlattice (or rather the square root of that)
Nxy = 64; Nz = 150; Nsuper = 1;
zMax = 6; zMin = -2;%units Å

disp("M * [0,1] = ")
disp(const.sheerMat*[0;1])

V = zeros(Nxy,Nxy,Nz);
X = zeros(Nxy,Nxy);
Y = zeros(Nxy,Nxy);
Xsuper = zeros(Nxy*Nsuper,Nxy*Nsuper);
Ysuper = zeros(Nxy*Nsuper,Nxy*Nsuper);
Z = linspace(zMin,zMax,Nz);

for i = 1:Nxy
    for j = 1:Nxy
        X(i,j) = (A1(1)*i+A2(1)*j)./Nxy;
        Y(i,j) = (A1(2)*i+A2(2)*j)./Nxy;
    end
end

for i = 1:Nxy*Nsuper
    for j = 1:Nxy*Nsuper
        Xsuper(i,j) = (A1(1)*i+A2(1)*j)./Nxy;
        Ysuper(i,j) = (A1(2)*i+A2(2)*j)./Nxy;
    end
end
  for k = 1:Nz
        V(:,:,k) = Vfunc(X,Y,Z(k));
  end
%We strictly ought to be careful with boundary conditions cos MS doesn't
%actually check them lol
%===
%% Now we duplicate the lattice to get a superlattice

Vsuper = zeros(Nsuper*Nxy,Nsuper*Nxy,Nz);
for z = 1:Nz
    for nx = 1:Nxy:Nsuper*Nxy
        for ny = 1:Nxy:Nsuper*Nxy
            Vsuper(nx:nx+Nxy-1,ny:ny+Nxy-1,z) = V(:,:,z);
        end
    end
end
%writematrix(Vsuper,"V.csv")

%Vsuper = readmatrix("V_boyao.csv");
%Vsuper = reshape(Vsuper,[Ncell,Ncell,Nz]);

%===
%% Now add imperfections to the lattice TODO update this for MoS2

%===
%% We also prepare a .csv which contains an equipotential plot.
equipotValue = 0;%Units meV ig
eqCharArr = [num2str(equipotValue,'%+g') , ' meV'];
equipotentialMat = zeros(Nxy*Nsuper,Nxy*Nsuper);
M = max(Vsuper,[],"all");
for i = 1:Nxy*Nsuper
    for j = 1:Nxy*Nsuper
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

%==

[S_X, S_Y, Mo_X, Mo_Y] = GetLatticePoints(Nsuper);

scatter(S_X,S_Y,'filled','o',SizeData=400)
hold on
daspect([1 1 1])
scatter(Mo_X,Mo_Y,'filled','pentagram',SizeData=400)
daspect([1 1 1])
hold off
%===
%% Plot the potential
% Plot of a slice of the potential in the nth row, that is for constant x
  row = floor(Nxy/2);
figure
contourf(Z,  linspace(0, sqrt(3)*const.b*Nsuper, Nxy*Nsuper), ...%!!!
    reshape(Vsuper(row,:,:), [Nxy*Nsuper,Nz]), linspace(-20,100,24))

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
zRow = floor((zSample - zMin)/(zMax-zMin) * Nz);
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

potStructArray.a1=A1; potStructArray.a2=A2;
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
        Q = cos(2*pi*x/const.a) + cos(2*pi*y/const.a);
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
        X_n = X ./ (const.b*3);
        Y_n = Y ./ (const.b*3/sqrt(3));
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
        %+ V1func(Z) * Qfunc(X,Y)...
    VmatrixElement = V0func(Z) ...
        + V1func(Z) * Qhexfunc(X,Y)...
        + Qhexfunc(X-const.b,Y) * V2func(Z);
end

function [DV] = Dropoff(z)
  %Use this to attenuate the gaussian in z
    DV = -exp(2*const.alpha*(const.z0-z));
end

function [X1coords,Y1coords,X2coords,Y2coords] = GetLatticePoints(Nsuper)
%Assumes 6 unit cells per lattice
  a1 = [const.b, 0];
  a2 = [const.b/2,const.b * sqrt(3)/2];
  X = zeros(2,3);
  Y = zeros(2,3);
  X1coords = zeros(Nsuper,3*Nsuper);
  Y1coords = zeros(Nsuper,3*Nsuper);
  X2coords = zeros(Nsuper,3*Nsuper);
  Y2coords = zeros(Nsuper,3*Nsuper);
      X(1,1) = 0;Y(1,1) = 0;
      X(2,1) = a1(1);Y(2,1) = 0;
      X(1,2) = a1(1)+a2(1);Y(1,2) = a1(2)+a2(2);
      X(2,2) = a1(1)*2+a2(1);Y(2,2) = a1(2)*2+a2(2);
      X(1,3) = a1(1)*2+a2(1)*2;Y(1,3)=a1(2)*2+a2(2)*2;
      X(2,3) = a2(1)*2;Y(2,3)=a2(2)*2;
    for i = 1:Nsuper
    for j = 1:Nsuper
      disp("(i-1)*Nsuper = " + (i-1)*Nsuper + ", (j-1)*Nsuper = " + (j-1)*Nsuper)
      bottomLeft = ((i-1)*a1+(j-1)*a2)*3;
      for l = 1:3
        X1coords((i-1)+1,(j-1)*3+l) = X(1,l) + bottomLeft(1);
        Y1coords((i-1)+1,(j-1)*3+l) = Y(1,l) + bottomLeft(2);
        X2coords((i-1)+1,(j-1)*3+l) = X(2,l) + bottomLeft(1);
        Y2coords((i-1)+1,(j-1)*3+l) = Y(2,l) + bottomLeft(2);
      end
    end
    end
end

function [Vout] = AddDefect(Vin,Nsuper)
  Vout = Vin;
  Nxy = size(Vout,1);
  Nz = size(Vout,3);
  for k = 1:Nz
    dropoff = Dropoff(Z(k));
    for i = 1:Nxy*Nsuper
      for j = 1:Nxy*Nsuper
        x = Xsuper(i,j);
        y = Ysuper(i,j);
        val = Gaussian2D(x,y, ...
          [const.a*1.3,const.a*1.3],const.a/2,3*const.D*dropoff);
        Vout(i,j,k) = Vout(i,j,k)+val;
        disp("x, y, z = " + x + ", " + y + ", " + Z(k) +...
            ", Value = " + val);
      end
    end
  end
end

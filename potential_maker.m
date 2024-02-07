clear; close all; clc;
rng default;

%Number of grid points, number of Z points, and number of lattices
%contained in the overall superlattice (or rather the square root of that)
Nxy = 64; Nz = 150; Nsuper = 1;
zMax = 6; zMin = -2;%units Å

%a = 2.84Å. see const.m for more stuff
%a1=[const.a,0];
%a2=[0,const.a]; 
a1=[const.c,0];
a2=[const.c/2,const.c * sqrt(3)/2];
a3=[0,0,const.c];
%A1 = a1;
%A2 = a2;
[b1,b2,b3] = Reciprocal([a1,0],[a2,0],a3);
%% data for python hex plotter WIP
writematrix([],'latticeVects.info_for_vivian_python_nice_plotting_hexagon_script',FileType='text')
a1str = [char(num2str(a1))];
a2str = [char(num2str(a2))];
b1str = [char(num2str(b1(1:2)))];
b2str = [char(num2str(b2(1:2)))];
S = fileread('latticeVects.info_for_vivian_python_nice_plotting_hexagon_script');
realStr = ['Real space vectors:',newline,'a1 = ',a1str, newline, 'a2 = ',a2str];
recpStr = ['Reciprocal vectors:',newline,'b1 = ',b1str, newline, 'b2 = ', b2str];

S = [realStr,newline,recpStr,S];
FID = fopen('latticeVects.info_for_vivian_python_nice_plotting_hexagon_script', 'w');
if FID == -1, error('Cannot open file %s', FileName); end
fwrite(FID, S, 'char');
fclose(FID);
%%



V = zeros(Nxy,Nxy,Nz);
X = zeros(Nxy,Nxy);
Y = zeros(Nxy,Nxy);
Xsuper = zeros(Nxy*Nsuper,Nxy*Nsuper);
Ysuper = zeros(Nxy*Nsuper,Nxy*Nsuper);
Z = linspace(zMin,zMax,Nz);

for i = 1:Nxy
    for j = 1:Nxy
        X(i,j) = (a1(1)*i+a2(1)*j)./Nxy;
        Y(i,j) = (a1(2)*i+a2(2)*j)./Nxy;
    end
end

for i = 1:Nxy*Nsuper
    for j = 1:Nxy*Nsuper
        Xsuper(i,j) = (a1(1)*i+a2(1)*j)./Nxy;
        Ysuper(i,j) = (a1(2)*i+a2(2)*j)./Nxy;
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
%% Now add imperfections to the lattice
%Vsuper = AddSulphurDefect(Vsuper,2,4,Xsuper,Ysuper,Z,Nsuper);
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
                %disp("i, j, k = " + i + ", " + j + ", " + k +...
                %    ", Value = " + Vsuper(i,j,k) +...
                %    ", Z(k) = " + Z(k));
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
contourf(Z,  linspace(0, const.c*Nsuper, Nxy*Nsuper), ...%!!!
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
    function [V0] = V0func(z,z0)
        V0 = const.D * exp(2*const.alpha*(z0-z))...
            -2*const.D*exp(const.alpha*(z0-z));
    end
    function [V1] = V1func(z,z0,wellDepth)
        V1 = 2*const.beta*wellDepth*exp(2*const.alpha*(z0-z));
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
        x_n = x / (const.c);
        y_n = y / (const.c/sqrt(3));
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
        X_n = X ./ (const.c);
        Y_n = Y ./ (const.c*sqrt(3));
        Q = zeros(length(X),length(Y));
        
        mu = X_n - Y_n;
        nu = Y_n*2;
        Q = Q + cos(2*pi*nu) + cos(2*pi*mu);

        mu = Y_n.*2;
        nu = -Y_n - X_n;
        Q = Q + cos(2*pi*nu) + cos(2*pi*mu);

        nu = Y_n - X_n;
        mu = -Y_n - X_n;
        Q = Q + cos(2*pi*nu) + cos(2*pi*mu);
        
        Q = Q./3;
        %Q = cos(2*pi*nu/const.a)^5 + cos(2*pi*mu/const.a)^5;
    end
        %+ V1func(Z) * Qfunc(X,Y)...
    VmatrixElement = V0func(Z,1) ...
        + V1func(Z,1.5,7.63) * Qhexfunc(X,Y)...
        + Qhexfunc(X-const.c/2,Y-(const.c*1/(2*sqrt(3)))) * V1func(Z,1,7.63);
      %VmatrixElement = Qhexfunc(X,Y) * Dropoff(Z) * const.D;
end

function [DV] = Dropoff(z)
  %Use this to attenuate the gaussian in z
    DV = exp(2*const.alpha*(const.z0-z));
end

function [X1coords,Y1coords,X2coords,Y2coords] = GetLatticePoints(Nsuper)
%Assumes 6 unit cells per lattice
  a1 = [const.c, 0];
  a2 = [const.c/2,const.c * sqrt(3)/2];
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

function [Vout] = AddSulphurDefect(Vin,m1,m2,Xsuper,Ysuper,Z,Nsuper)
%Adds a defect at sulphur site (m1,m2)
  Vout = Vin;
  NxySuper = size(Vout,1);
  Nz = size(Vout,3);
  [X1,Y1,X2,Y2] = GetLatticePoints(Nsuper);
  centre = [X1(m1,m2),Y1(m1,m2)];
  disp("Centre:")
  disp(centre)
  for k = 1:Nz
    dropoff = Dropoff(Z(k));
    for i = 1:NxySuper
      for j = 1:NxySuper
        x = Xsuper(i,j);
        y = Ysuper(i,j);
        val = Gaussian2D(x,y, ...
          centre,const.c/2,const.D*dropoff*1.1);
        %I made these values the fuck up
        Vout(i,j,k) = Vout(i,j,k)+val;
        %disp("x, y, z = " + x + ", " + y + ", " + Z(k) +...
        %    ", Value = " + val);
      end
    end
  end
end

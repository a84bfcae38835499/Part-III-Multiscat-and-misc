clear; close all; clc;
rng default;

%Number of grid points, number of Z points, and number of lattices
%contained in the overall superlattice (or rather the square root of that)
Nxy = 64; Nz = 150; Nsuper = 2;
zMax = 8; zMin = -1;%units Å

%a = 2.84Å. see const.m for more stuff
%a1=[const.a,0];
%a2=[0,const.a]; 
%a1=[const.c,0];
%a2=[const.c/2,const.c * sqrt(3)/2];
a1=[-const.c,0];
a2=[const.c/2,const.c*sqrt(3)/2];
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
%Vsuper = AddSulphurDefect(Vsuper,2,1,a1,a2,Xsuper,Ysuper,Z);
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
    shading interp

%% Plot the potential
zSample = 2.0;
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

potStructArray.a1=Nsuper*a1; potStructArray.a2=Nsuper*a2;
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
  function [V0] = V0func(z,z0,backgroundDepth)
        V0 = backgroundDepth * exp(2*const.alpha*(z0-z))...
            -2*backgroundDepth*exp(const.alpha*(z0-z));
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
        Q = 2/3 * (cos(2*pi*(X_n-Y_n))+cos(4*pi*Y_n)+cos(2*pi*(X_n+Y_n)));
        %Q = cos(2*pi*nu/const.a)^5 + cos(2*pi*mu/const.a)^5;
    end
        %+ V1func(Z) * Qfunc(X,Y)...
    VmatrixElement = V0func(Z,2.4,const.MoS2Depth/3) ...
        + V1func(Z,3,const.MoS2Depth/3) * Qhexfunc(X,Y)...
        + Qhexfunc(X-const.c/2,Y-(const.c*1/(2*sqrt(3)))) * V1func(Z,1.7,const.MoS2Depth/3);
      %VmatrixElement = Qhexfunc(X,Y) * Dropoff(Z) * const.D;
end

function [DV] = Dropoff(z,z0)
  %Use this to attenuate the gaussian in z
    DV = exp(2*const.alpha*(z0-z));
end

function [Vout] = AddSulphurDefect(Vin,m1,m2,a1,a2,Xsuper,Ysuper,Z)
%Adds a defect at sulphur site (m1,m2)
  Vout = Vin;
  NxySuper = size(Vout,1);
  Nz = size(Vout,3);
  centre = m1*a1+m2*a2;
  disp("Centre:")
  disp(centre)
  for k = 1:Nz
    for i = 1:NxySuper
      for j = 1:NxySuper
        x = Xsuper(i,j);
        y = Ysuper(i,j);
        val = Gaussian2D(x,y, ...
          centre,const.c*0.2,-3*const.beta*const.MoS2Depth*exp(2*const.alpha*(3-Z(k))));
        %These parameters have been fined tuned to match the requirements
        Vout(i,j,k) = Vout(i,j,k)+val;
        %disp("x, y, z = " + x + ", " + y + ", " + Z(k) +...
        %    ", Value = " + val);
      end
    end
  end
end
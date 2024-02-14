clear; close all; clc;
rng default;

%Number of grid points, number of Z points, and number of lattices
%contained in the overall superlattice (or rather the square root of that)
Nxy = 64; Nz = 100; Nsuper = 2;
zMax = 6; zMin = -1;%units Å

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

for i = 0:Nxy-1
    for j = 0:Nxy-1
        X(i+1,j+1) = (a1(1)*i+a2(1)*j)./Nxy;
        Y(i+1,j+1) = (a1(2)*i+a2(2)*j)./Nxy;
    end
end

for i = 0:Nxy*Nsuper-1
    for j = 0:Nxy*Nsuper-1
        Xsuper(i+1,j+1) = (a1(1)*i+a2(1)*j)./Nxy;
        Ysuper(i+1,j+1) = (a1(2)*i+a2(2)*j)./Nxy;
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
%Vsuper = AddSulphurDefect(false,Vsuper,1,1,a1,a2,Nsuper,Xsuper,Ysuper,Z);
%===
%% Plot the potential
%nPlot = 2/3;mPlot = 1/2;
nPlotDef = 1;mPlotDef = 1;
aboveCol = [0.8 0.3 1];
PlotPotentialAlongZ(Vsuper,a1,a2,mPlotDef,nPlotDef,Z,0,dft.aboveDefect,aboveCol)
nPlotHol = 2/3;mPlotHol = 1/3;
holCol = [0.0 1 0.6];
PlotPotentialAlongZ(Vsuper,a1,a2,mPlotHol,nPlotHol,Z,0,dft.aboveHollow,holCol)
nPlotMo = 1/3;mPlotMo = 2/3;
moCol = [1 0.4 0];
PlotPotentialAlongZ(Vsuper,a1,a2,mPlotMo,nPlotMo,Z,0,dft.aboveMo,moCol)
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
xlim([2,6])
title('Potential in z, used in simulation')
hbar = colorbar;
ylabel(hbar,'Energy / meV');

%% Plot the potential
% Linearly interpolated equipotential plot
fontsize(gcf,scale=1)
%equipotential_plot('V', Vsuper, 'z', Z, 'X', Xsuper, 'Y', Ysuper)
shading interp
%hold on

importfile("DFT_Pure.mat")
x1=[const.d,0];
x2=[const.d/2,const.d * sqrt(3)/2];
DFTSuper = 1;
XDFTSuper = zeros(12*DFTSuper);
YDFTSuper = zeros(12*DFTSuper);

for i = 0:12*DFTSuper-1
    for j = 0:12*DFTSuper-1
        XDFTSuper(i+1,j+1) = (x1(1)*i+x2(1)*j)./12;
        YDFTSuper(i+1,j+1) = (x1(2)*i+x2(2)*j)./12;
    end
end


VDFTsuper = zeros(DFTSuper*12,DFTSuper*12,19);
for z = 1:19
    for nx = 1:12:DFTSuper*12
        for ny = 1:12:DFTSuper*12
            VDFTsuper(nx:nx+12-1,ny:ny+12-1,z) = pagetranspose(Pot_M(:,:,z));
        end
    end
end

equipotential_plot('V',VDFTsuper,'z',Z(1:19),'X',XDFTSuper,'Y',YDFTSuper)
%hold off
%% Plot the potential
zSample = 2.0;
zRow = floor((zSample - zMin)/(zMax-zMin) * Nz);
figure
contourf(Xsuper,Ysuper,Vsuper(:,:,zRow),10)
daspect([1 1 1])
xlabel('x/Å')
ylabel('y/Å')
title('Potentials at z = ' + string(zSample) + ' Å');
colormap(parula(15))
hbar = colorbar('southoutside');
xlabel(hbar,'Energy / meV');

    fontsize(gcf,scale=1)
hold on
    xPlot = mPlotDef*a1(1)+nPlotDef*a2(1);
    yPlot = mPlotDef*a1(2)+nPlotDef*a2(2);
    plot(xPlot,yPlot,'*',MarkerSize=24,Color=aboveCol);
    plot(xPlot,yPlot,'.',MarkerSize=24,Color=aboveCol);

    xPlot = mPlotHol*a1(1)+nPlotHol*a2(1);
    yPlot = mPlotHol*a1(2)+nPlotHol*a2(2);
    plot(xPlot,yPlot,'*',MarkerSize=24,Color=holCol);
    plot(xPlot,yPlot,'.',MarkerSize=24,Color=holCol);

    xPlot = mPlotMo*a1(1)+nPlotMo*a2(1);
    yPlot = mPlotMo*a1(2)+nPlotMo*a2(2);
    plot(xPlot,yPlot,'*',MarkerSize=24,Color=moCol);
    plot(xPlot,yPlot,'.',MarkerSize=24,Color=moCol);
 hold off
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
    VmatrixElement = V0func(Z,const.zOffset+2.4,const.MoS2Depth/3) ...
        + V1func(Z,const.zOffset+3,const.MoS2Depth/3) * Qhexfunc(X,Y)...
        + Qhexfunc(X-const.c/2,Y-(const.c*1/(2*sqrt(3)))) * V1func(Z,const.zOffset+1.7,const.MoS2Depth/3);
      %VmatrixElement = Qhexfunc(X,Y) * Dropoff(Z) * const.D;
end

function [DV] = Dropoff(z,z0)
  %Use this to attenuate the gaussian in z
    DV = exp(2*const.alpha*(z0-z));
end

function [Vout] = AddSulphurDefect(doWeRepeat,Vin,m,n,a1,a2,Nsuper,Xsuper,Ysuper,Z)
%Adds a defect at sulphur site (m1,m2)
  Vout = Vin;
  NxySuper = size(Vout,1);
  Nz = size(Vout,3);
  centresX = zeros(3);
  centresY = zeros(3);
  centre = m*a1+n*a2;
  disp("Centre:")
  disp(centre)
  if doWeRepeat %This is awful optimisation lolololo
    for m = -1:1
      for n = -1:1
        centresX(m+2,n+2) = centre(1)+m*a1(1)*Nsuper+n*a2(1)*Nsuper;
        centresY(m+2,n+2) = centre(2)+m*a1(2)*Nsuper+n*a2(2)*Nsuper;
      end
    end
  end
  if doWeRepeat
    for m = -1:1
      for n = -1:1
        for k = 1:Nz
          for i = 1:NxySuper
            for j = 1:NxySuper
              x = Xsuper(i,j);
              y = Ysuper(i,j);
              centre = [centresX(m+2,n+2) centresY(m+2,n+2)];
              val = Gaussian2D(x,y, ...
                centre,const.c*0.2,-3*const.beta*const.MoS2Depth* ...
                exp(2*const.alpha*(const.zOffset+3-Z(k))));
              %These parameters have been fined tuned to match the requirements
              Vout(i,j,k) = Vout(i,j,k)+val;
              %disp("x, y, z = " + x + ", " + y + ", " + Z(k) +...
              %    ", Value = " + val);
            end
          end
        end
      end
    end
  else
    for k = 1:Nz
      for i = 1:NxySuper
        for j = 1:NxySuper
          x = Xsuper(i,j);
          y = Ysuper(i,j);
          val = Gaussian2D(x,y, ...
            centre,const.c*0.2,-3*const.beta*const.MoS2Depth* ...
            exp(2*const.alpha*(3+const.zOffset-Z(k))));
          %These parameters have been fined tuned to match the requirements
          Vout(i,j,k) = Vout(i,j,k)+val;
          %disp("x, y, z = " + x + ", " + y + ", " + Z(k) +...
          %    ", Value = " + val);
        end
      end
    end
  end
end

function PlotPotentialAlongZ(V,a1,a2,m,n,Z,zMin,dftPot,plotColor)
  NxyNsuper = size(V,1);
  SpaghettiBolognaise = [a1(1) a2(1);a1(2) a2(2)]/NxyNsuper;
  k = int8(interp1(Z,1:numel(Z),zMin));

  centre = m*a1+n*a2;
  result = SpaghettiBolognaise\(centre');
  i = int8(result(1));
  j = int8(result(2));
  Vpiece = squeeze(V(i,j,:));
  figure
  p = plot(Z(k:end),Vpiece(k:end),'DisplayName','Estimated analytical potential');
  p.LineStyle = ":";
  p.Color = plotColor;
  p.Marker = ".";
  tit = 'Plot of potential at x, y = ' + string(centre(1)) + ', ' + string(centre(2));
  title(tit)
  xlabel('z/Å')
  ylabel('Energy/meV')
  hold on
  legend()
  d = plot(dft.zAxis,dftPot,'DisplayName','DFT result');
  d.LineStyle = "-";
  d.Color = [0 0.4470 0.7410];
  d.Marker = ".";
  ylim([1.1*min(dftPot) 1.1*max(dftPot)])
  hold off
end
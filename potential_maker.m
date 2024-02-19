clear; close all; clc;
rng default;

%Number of grid points, number of Z points, and number of lattices
%contained in the overall superlattice (or rather the square root of that)
Nxy = 64; Nz = 50; Nsuper = 1;
zMax = 6; zMin = 1.5;%units Å

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
%% Import Min's DFT

importfile("DFT_Pure.mat")
x1=[const.d,0];
x2=[const.d/2,const.d * sqrt(3)/2];
DFTsuper = 2;
XDFTsuper = zeros(12*DFTsuper);
YDFTsuper = zeros(12*DFTsuper);
ZDFT = linspace(1.5,6,19);

for i = 0:12*DFTsuper-1
    for j = 0:12*DFTsuper-1
        XDFTsuper(i+1,j+1) = (x1(1)*i+x2(1)*j)./12;
        YDFTsuper(i+1,j+1) = (x1(2)*i+x2(2)*j)./12;
    end
end

XDFTsuper = XDFTsuper - const.c/(2-0.3); %makes the 0,0 point be a sulphur
XDFTsuper = XDFTsuper - x1(1) - x2(1);%shifts potential to overlap with the smaller unit cell, for interpolation purposes
YDFTsuper = YDFTsuper - x1(2) - x2(2);


theta = 30;
rotMat = [cosd(theta) -sind(theta);
          sind(theta)  cosd(theta)];

for i = 1:size(XDFTsuper,1)
  for j = 1:size(YDFTsuper,1)
    vIn = [XDFTsuper(i,j); YDFTsuper(i,j)];
    vOut = rotMat * vIn;
    XDFTsuper(i,j) = vOut(1);
    YDFTsuper(i,j) = vOut(2);
  end
end

VDFTsuper = zeros(DFTsuper*12,DFTsuper*12,19);
for z = 1:19
    for nx = 1:12:DFTsuper*12
        for ny = 1:12:DFTsuper*12
            VDFTsuper(nx:nx+12-1,ny:ny+12-1,z) = pagetranspose(Pot_M(:,:,z))*1000;
        end
    end
end

[y1,y2,y3] = Reciprocal([x1,0],[x2,0],a3);

%%
V = zeros(Nxy,Nxy,Nz);
Vinterp = zeros(Nxy,Nxy,Nz);
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
%% Now interpolate the DFT data into a useful basis
interpolateDFTdata = true;
Vvect = zeros(Nz*Nxy*Nxy,1);
if interpolateDFTdata == true
  VDFTvect = zeros(DFTsuper*DFTsuper*12*12*19,1);
  XDFTvect = VDFTvect;
  YDFTvect = VDFTvect;
  ZDFTvect = VDFTvect;
  index = 0;
  for k = 1:19
    z = ZDFT(k);
     for j = 1:12*DFTsuper
      for i = 1:12*DFTsuper
        if(index + 1 ~= 144*DFTsuper*DFTsuper*(k-1)+12*DFTsuper*(j-1)+i)
          error("F")
        end
        index = 144*DFTsuper*DFTsuper*(k-1)+12*DFTsuper*(j-1)+i;
        %disp("index = " + num2str(index))
        XDFTvect(index) = XDFTsuper(i,j);
        YDFTvect(index) = YDFTsuper(i,j);
        ZDFTvect(index) = z;
        VDFTvect(index) = VDFTsuper(i,j,k);
      end
    end
  end
  InterpolatedFn = scatteredInterpolant(XDFTvect,YDFTvect,ZDFTvect,VDFTvect,'natural','none');

  Xvect = squeeze(zeros(Nxy*Nxy*Nz,1));
  Yvect = Xvect;
  Zvect = Xvect;
  for k = 1:Nz
    z = Z(k);
     for j = 1:Nxy
      for i = 1:Nxy
        index2 = Nxy*Nxy*(k-1)+Nxy*(j-1)+i;
        %disp("index2 = " + num2str(index2))
        Xvect(index2) = X(i,j);
        Yvect(index2) = Y(i,j);
        Zvect(index2) = z;
      end
    end
  end

  Vvect = InterpolatedFn(Xvect,Yvect,Zvect);
  if(anynan(Vvect))
    error("Nan found!")
  end

  for k = 1:Nz
    for j = 1:Nxy
      for i = 1:Nxy
        V(i,j,k) = Vvect(Nxy*Nxy*(k-1)+Nxy*(j-1)+i);
      end   
    end
  end
end
%% Now we duplicate the lattice to get a superlattice

Vsuper = zeros(Nsuper*Nxy,Nsuper*Nxy,Nz);
for z = 1:Nz
    for nx = 1:Nxy:Nsuper*Nxy
        for ny = 1:Nxy:Nsuper*Nxy
            Vsuper(nx:nx+Nxy-1,ny:ny+Nxy-1,z) = V(:,:,z);
        end
    end
end
Vinterpsuper = zeros(Nsuper*Nxy,Nsuper*Nxy,Nz);
for z = 1:Nz
    for nx = 1:Nxy:Nsuper*Nxy
        for ny = 1:Nxy:Nsuper*Nxy
            Vinterpsuper(nx:nx+Nxy-1,ny:ny+Nxy-1,z) = Vinterp(:,:,z);
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
%% Plot the potential. Disabled for now, as if the grid res is too high it complains
%nPlot = 2/3;mPlot = 1/2;
nPlotDef = 1;mPlotDef = 1;
aboveCol = [0.3 0. 1];
ComparePotentials(Vsuper,Vinterpsuper,'Analytical potential','DFT interpolated',a1,a2,mPlotDef,nPlotDef,Z,0,aboveCol)
nPlotHol = 2/3;mPlotHol = 1/3;
holCol = [0.0 0.6 0.2];
ComparePotentials(Vsuper,Vinterpsuper,'Analytical potential','DFT interpolated',a1,a2,mPlotHol,nPlotHol,Z,0,holCol)
nPlotMo = 1/3;mPlotMo = 2/3;
moCol = [1 0.2 0];
ComparePotentials(Vsuper,Vinterpsuper,'Analytical potential','DFT interpolated',a1,a2,mPlotMo,nPlotMo,Z,0,moCol)

%% Get min and max bounds of the potentials
DFTmin = min(VDFTsuper,[],"all")
DFTmax = max(VDFTsuper,[],"all")
AnalyticMin = min(Vsuper,[],"all")
AnalyticMax = max(Vsuper,[],"all")

%% Now change all the crap to be Min's DFT
copyDFT = false;
if copyDFT
  Nsuper = DFTsuper;
  Nxy = 12;
  Nz = 19;
  Xsuper = XDFTsuper;
  Ysuper = YDFTsuper;
  Z = ZDFT;
  Vsuper = VDFTsuper;
  a1 = x1;a2=x2;
  b1 = y1; b2 = y2;
end
%% data for python hex plotter
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

%% Plot the potential
% Plot of a slice of the potential in the nth row, that is for constant x
  row = floor(Nxy/2);
figure
contourf(Z,  linspace(0, const.c*Nsuper, Nxy*Nsuper), ...%!!!
    reshape(Vsuper(row,:,:), [Nxy*Nsuper,Nz]), linspace(-30,100,24))

    fontsize(gcf,scale=1)
xlabel('z/Å')
ylabel('y/Å') %is this x or y? I think y but idrk
colorbar
xlim([1.5,6])
title('Potential in z, used in simulation')
hbar = colorbar;
ylabel(hbar,'Energy / meV');
figure
fileindx = 1;
for i = 0
  Vsoup = single(i);
  equipotential_plot('V', Vsuper, 'V0', Vsoup, 'z', Z, 'X', Xsuper, 'Y', Ysuper)
  shading interp
  hold on
  view([15 45])
  %equipotential_plot('V',VDFTsuper,'V0', Vsoup, 'z',ZDFT,'X',XDFTsuper,'Y',YDFTsuper)
  shading interp
  xlim([-3.5 2]*Nsuper);
  ylim([-0.5 3]*Nsuper);
  daspect([1 1 1])
  hold off
  savestr = "Figures/Frames/frame_" +num2str(fileindx,'%06d')+ ".jpg"
  fileindx = fileindx + 1;
  saveas(gcf,savestr,'jpg')
  figure
  %clf
end
%% Plot the potential
fontsize(gcf,scale=1)
zSample = 3;
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
%add indicators for where we're sampling the potential z
fontsize(gcf,scale=1)
plotPoints = false;
if(plotPoints)
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
end
%===

%% We supply the lattice to the mulitscat script so it can do its thing
doingMSshit = true;
if(doingMSshit)
    potStructArray.V = Vsuper;
    Multiscat.PreparePotentialFiles(potStructArray);
    
    Multiscat.prepareFourierLabels(Vsuper);
    
    potStructArray.a1=Nsuper*a1; potStructArray.a2=Nsuper*a2;
    potStructArray.zmin=Z(1);
    potStructArray.zmax=Z(end);
    potStructArray.zPoints=length(Z);
    
    confStruct=Multiscat.createConfigStruct(potStructArray);
    Multiscat.prepareConfigFile(confStruct);
end
%===

%[X,Y] = meshgrid(-1:0.01:1,-1:0.01:1);
%Y_prime = Y/sqrt(3);
%Z = ((cos(2*pi*(X-Y_prime))+cos(4*pi*Y_prime)+cos(2*pi*(X+Y_prime))) + 3/2)/(4.5);
%surf(X,Y,Z)
%shading('interp')
%colormap(plasma)

%% Function definitions

function [b1,b2,b3] = Reciprocal(a1,a2,a3)
    factor = 2*pi/dot(cross(a1,a2),a3);
    b1 = factor*cross(a2,a3);
    b2 = factor*cross(a3,a1);
    b3 = factor*cross(a1,a2);
end

function [VmatrixElement] = Vfunc(X,Y,Z)
  function [V0] = V0func(z,z0,backgroundDepth,alpha)
        V0 = backgroundDepth * exp(2*alpha*(z0-z))...
            -2*backgroundDepth*exp(alpha*(z0-z));
    end
    function [V1] = V1func(z,z0,D,alpha)
        V1 = 2*const.beta*D*exp(2*alpha*(z0-z));
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
        Q = ((cos(2*pi*(X_n-Y_n))+cos(4*pi*Y_n)+cos(2*pi*(X_n+Y_n))) + 3/2)/(4.5);
        %Q = cos(2*pi*nu/const.a)^5 + cos(2*pi*mu/const.a)^5;
    end
        %+ V1func(Z) * Qfunc(X,Y)...
    VmatrixElement = (V0func(Z,const.zOffset+2.1,25,1.2) ...
       + V1func(Z,const.zOffset+3.7,0.5,1.6))... %blue
       * Qhexfunc(X,Y) ...
       + (V0func(Z,const.zOffset+2.15,20,1.2) + ...
      + V1func(Z,const.zOffset+3,0,1.1)) ... % green
      * Qhexfunc(X-const.c/2,Y-(const.c*1/(2*sqrt(3)))) ...
      + (V0func(Z,const.zOffset+2.1,23,1.2) ...
      + V1func(Z,const.zOffset+1,15,1.1)) ... %red
      * Qhexfunc(X,Y - (const.c/sqrt(3)));
      %VmatrixElement = Qhexfunc(X,Y) * Dropoff(Z) * const.D;
end

function [DV] = Dropoff(z,z0)
  %Use this to attenuate the gaussian in z
    DV = exp(2*const.alpha*(z0-z));
end

function [Vout] = AddSulphurDefect(doWeRepeat,Vin,m,n,a1,a2,Nsuper,Xsuper,Ysuper,Z)
factor = -0.7*const.beta*const.MoS2Depth;
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
                centre,const.c*0.2, ...
                factor * exp(2*const.alpha*(const.zOffset+3-Z(k))));
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
            centre,const.c*0.2, ...
            factor * exp(2*const.alpha*(3+const.zOffset-Z(k))));
          %These parameters have been fined tuned to match the requirements
          Vout(i,j,k) = Vout(i,j,k)+val;
          %disp("x, y, z = " + x + ", " + y + ", " + Z(k) +...
          %    ", Value = " + val);
        end
      end
    end
  end
end

function ComparePotentials(V1,V2,V1name,V2name,a1,a2,m,n,Z,zMin,plotColor)
  NxyNsuper = size(V1,1);
  SpaghettiBolognaise = [a1(1) a2(1);a1(2) a2(2)]/NxyNsuper;
  k = int8(interp1(Z,1:numel(Z),zMin));
  if(k == 0)
    k = 1;
  end
  centre = m*a1+n*a2;
  result = SpaghettiBolognaise\(centre');
  i = int8(result(1));
  j = int8(result(2));
  V1piece = squeeze(V1(i,j,:));
  V2piece = squeeze(V2(i,j,:));
  disp([i j k])
  figure
  p = plot(Z(k:end),V1piece(k:end),'DisplayName',V1name);
  p.LineStyle = ":";
  p.Color = plotColor;
  p.Marker = ".";
  p.LineWidth=2;
  tit = 'Plot of potentials at x, y = ' + string(centre(1)) + ', ' + string(centre(2));
  title(tit)
  xlabel('z/Å')
  ylabel('Energy/meV')
  hold on
  legend()
  d = plot(Z(k:end),V2piece(k:end),'DisplayName',V2name);
  d.LineStyle = "-";
  d.Color = [0 0.4470 0.7410];
  d.Marker = ".";
  d.LineWidth=2;
  xlim([1.5 6])
  ylim([-30 100])
  hold off
end
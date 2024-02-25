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

Theta = 0.5;
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
        Vinterp(i,j,k) = Vvect(Nxy*Nxy*(k-1)+Nxy*(j-1)+i);
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
%% oh fuck
copyInterp = false;
if(copyInterp)
  Vsuper = Vinterpsuper;
end

%% Now add imperfections to the lattice
Nsites = Nsuper*Nsuper;
disp("Target number of sites = " + (Nsites * Theta))
Ndefect = int8(Nsites * Theta);
ms = 0:1:Nsites;
ns = 0:1:Nsites;

%The question is - How do we decide which sites to add defects on to?
%Do we directly implement avoidance of nearest-neighbors?
avoidNearestNeighbors = false
if(avoidNearestNeighbors)
  error("Nearest neighbour avoidance not yet implemented!")
else

end
%Vsuper = AddSulphurDefect(false,Vsuper,1,1,a1,a2,Nsuper,Xsuper,Ysuper,Z);
%===
%% Plot the potential. Disabled for now, as if the grid res is too high it complains
%nPlot = 2/3;mPlot = 1/2;
nPlotDef = 1;mPlotDef = 1;
aboveCol = [0.3 0. 1];
ComparePotentials(Vsuper,Vinterp,'Analytical potential','DFT interpolated',a1,a2,mPlotDef,nPlotDef,Z,Z,0,aboveCol)
nPlotHol = 2/3;mPlotHol = 1/3;
holCol = [0.0 0.6 0.2];
ComparePotentials(Vsuper,Vinterp,'Analytical potential','DFT interpolated',a1,a2,mPlotHol,nPlotHol,Z,Z,0,holCol)
nPlotMo = 1/3;mPlotMo = 2/3;
moCol = [1 0.2 0];
ComparePotentials(Vsuper,Vinterp,'Analytical potential','DFT interpolated',a1,a2,mPlotMo,nPlotMo,Z,Z,0,moCol)

%% for fitit
Zdefect = dft.zAxis;
Vdefect = dft.aboveDefect;
SpaghettiBolognaise = [a1(1) a2(1);a1(2) a2(2)]/(Nxy*Nsuper);
zFitMin = 1.5;
k = int8(interp1(Z,1:numel(Z),zFitMin));
if(k == 0)
  k = 1;
end
m = 1; n = 1;
centre = m*a1+n*a2;
result = SpaghettiBolognaise\(centre');
i = int8(result(1));
j = int8(result(2));
V1piece = squeeze(Vinterpsuper(i,j,k:end));
Zpiece = Z(k:end);

weights = Zpiece;
for k = 1:numel(Zpiece)
  weights(k) = 1/(exp((3-Zpiece(k))*7)+1);
end
figure
plot(Zpiece,weights)
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
nsupstr = [char(num2str(Nsuper))];
S = fileread('latticeVects.info_for_vivian_python_nice_plotting_hexagon_script');
realStr = ['Real space vectors:',newline,'a1 = ',a1str, newline, 'a2 = ',a2str,newline,'Nsuper = ',nsupstr];
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
  figure
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
plotPoints = true;
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
  function [V] = V1(z)
    D1 = 17;
    alpha1 = 1.4267;
    b1 = -0.1991;
    beta1 = 0.3116;
    c1 = 14.1590;
    lambda1 = -0.6126;
    z01 = 3.4411;
    z11 = 5.5355;
    z21 = 3.5000;
    V = D1*(exp(2*alpha1*(z01-z))-2*exp(alpha1*(z01-z))-2*(b1*10)*exp(2*beta1*(z11-z))-c1*exp(lambda1*-(z21-z)));
  end

  function [V] = V2(z)
    D2 = 21;
    alpha2 = 1.2000;
    b2 = 2.0000;
    beta2 = 1.3000;
    c2 = 0.0400;
    lambda2 = 0.6000;
    z02 = 3.1500;
    z12 = 0;
    z22 = 6.0000;
    V = D2*(exp(2*alpha2*(z02-z))-2*exp(alpha2*(z02-z))-2*(b2*10)*exp(2*beta2*(z12-z))-c2*exp(lambda2*-(z22-z)));
  end

  function [V] = V3(z)
    D3 = 20;
    alpha3 = 1.1706;
    b3 = 0.6386;
    beta3 = 0.7456;
    c3 = 0.0064;
    lambda3 = 0.7557;
    z03 = 3.1578;
    z13 = -0.0261;
    z23 = 3.5021;
    V = D3*(exp(2*alpha3*(z03-z))-2*exp(alpha3*(z03-z))-2*(b3*10)*exp(2*beta3*(z13-z))-c3*exp(lambda3*-(z23-z)));
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
fittingDFT = false;
if(fittingDFT)
        %+ V1func(Z) * Qfunc(X,Y)...
    VmatrixElement = V1(Z) ... %blue, sulphur
       * Qhexfunc(X,Y) ...
       + V2(Z) ... %red, molybdenum
      * Qhexfunc(X-const.c/2,Y-(const.c*1/(2*sqrt(3)))) ...
      + V3(Z) ...%green, hollow site
      * Qhexfunc(X,Y - (const.c/sqrt(3)));
else
    VmatrixElement = (V0func(Z,const.zOffset+2.1,25,1.2) ...
       + V1func(Z,const.zOffset+3.7,0.5,1.6))... %blue
       * Qhexfunc(X,Y) ...
       + (V0func(Z,const.zOffset+2.15,20,1.2) + ...
      + V1func(Z,const.zOffset+3,0,1.1)) ... % green
      * Qhexfunc(X-const.c/2,Y-(const.c*1/(2*sqrt(3)))) ...
      + (V0func(Z,const.zOffset+2.1,23,1.2) ...
      + V1func(Z,const.zOffset+1,15,1.1)) ... %red
      * Qhexfunc(X,Y - (const.c/sqrt(3)));
end
      %VmatrixElement = Qhexfunc(X,Y) * Dropoff(Z) * const.D;
end

function [DV] = Dropoff(z,z0)
  %Use this to attenuate the gaussian in z
    DV = exp(2*const.alpha*(z0-z));
end

function [Vout] = AddSulphurDefect(doWeRepeat,Vin,m,n,a1,a2,Nsuper,Xsuper,Ysuper,Z)
%Adds a defect at sulphur site (m1,m2)
  Vout = Vin;
  function [v] = val(x,y,Z,k,centre)
    d = 8.4;
    gamma = 1.1;
    v = -d * exp(2*gamma*(4-Z(k))) * ...
    Gaussian2D(x,y, ...
    centre,const.c*0.2);
  end
  NxySuper = size(Vout,1);
  Nz = size(Vout,3);
  centresX = zeros(3);
  centresY = zeros(3);
  centre = m*a1+n*a2;
  disp("Centre:")
  disp(centre)

  if doWeRepeat
    for m = -1:1
      for n = -1:1
        centresX(m+2,n+2) = centre(1)+m*a1(1)*Nsuper+n*a2(1)*Nsuper;
        centresY(m+2,n+2) = centre(2)+m*a1(2)*Nsuper+n*a2(2)*Nsuper;
        for k = 1:Nz
          for i = 1:NxySuper
            for j = 1:NxySuper
              x = Xsuper(i,j);
              y = Ysuper(i,j);
              centre = [centresX(m+2,n+2) centresY(m+2,n+2)];
              Vout(i,j,k) = Vout(i,j,k)+val(x,y,Z,k,centre);
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
          Vout(i,j,k) = Vout(i,j,k)+val(x,y,Z,k,centre);
          %disp("x, y, z = " + x + ", " + y + ", " + Z(k) +...
          %    ", Value = " + val);
        end
      end
    end
  end
end

function ComparePotentials(V1,V2,V1name,V2name,a1,a2,m,n,Z1,Z2,zMin,plotColor)
  NxyNsuper1 = size(V1,1);
  SpaghettiBolognaise1 = [a1(1) a2(1);a1(2) a2(2)]/NxyNsuper1;
  k1 = int8(interp1(Z1,1:numel(Z1),zMin));
  if(k1 == 0)
    k1 = 1;
  end
  centre = m*a1+n*a2;
  result = SpaghettiBolognaise1\(centre');
  i1 = int8(result(1));
  j1 = int8(result(2));
  V1piece = squeeze(V1(i1,j1,:));
  
  NxyNsuper2 = size(V2,1);
  SpaghettiBolognaise2 = [a1(1) a2(1);a1(2) a2(2)]/NxyNsuper2;
  k2 = int8(interp1(Z2,1:numel(Z2),zMin));
  if(k2 == 0)
    k2 = 1;
  end
  centre = m*a1+n*a2;
  result = SpaghettiBolognaise2\(centre');
  i2 = int8(result(1));
  j2 = int8(result(2));
  V2piece = squeeze(V2(i2,j2,:));
  %V2piece = V2;
  disp([i1 j1 k1])
  disp([i2 j2 k2])
  figure
  p = plot(Z1(k1:end),V1piece(k1:end),'DisplayName',V1name);
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
  disp(size(Z2))
  disp(size(V2piece))
  d = plot(Z2(k2:end),V2piece(k2:end),'DisplayName',V2name);
  d.LineStyle = "-";
  d.Color = [0 0.4470 0.7410];
  d.Marker = ".";
  d.LineWidth=2;
  xlim([1.5 6])
  ylim([-30 100])
  hold off
end
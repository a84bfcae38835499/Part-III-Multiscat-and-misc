clear; close all; clc;
rng default;

%Number of grid points, number of Z points, and number of lattices
%contained in the overall superlattice (or rather the square root of that)
Nxy = 32; Nz = 100; Nsuper = 1;
Theta = 0.;
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
N1DFT = linspace(0,const.c,12*DFTsuper);
N2DFT = linspace(0,const.c,12*DFTsuper);
N1DFT_a = [N1DFT const.c*(12*DFTsuper+1/(12*DFTsuper))];
N2DFT_a = [N2DFT const.c*(12*DFTsuper+1/(12*DFTsuper))];
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

VDFT = zeros(12,12,19);
for k = 1:19
  VDFT(:,:,k) = pagetranspose(Pot_M(:,:,k))*1000;
end
VDFT_a = zeros(12+1,12+1,19);
VDFT_a(1:12,1:12,:) = VDFT;
VDFT_a(12+1,1:12,:) = VDFT(1,:,:);
VDFT_a(1:12,12+1,:) = VDFT(:,1,:);
VDFT_a(end ,end ,:) = VDFT(end,end,:);

VDFTsuper = zeros(DFTsuper*12,DFTsuper*12,19);
for z = 1:19
    for nx = 1:12:DFTsuper*12
        for ny = 1:12:DFTsuper*12
            VDFTsuper(nx:nx+12-1,ny:ny+12-1,z) = pagetranspose(Pot_M(:,:,z))*1000;
        end
    end
end

VDFTsuper_a = zeros(DFTsuper*12+1,DFTsuper*12+1,19);
VDFTsuper_a(1:DFTsuper*12,1:DFTsuper*12  ,:) = VDFTsuper;
VDFTsuper_a(DFTsuper*12+1,1:DFTsuper*12  ,:) = VDFTsuper(1,:,:);
VDFTsuper_a(1:DFTsuper*12,DFTsuper*12+1  ,:) = VDFTsuper(:,1,:);
VDFTsuper_a(end          ,end            ,:) = VDFTsuper(end,end,:);

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

if(interpolateDFTdata)
  oldmethod = true;
  if(oldmethod)
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
    %Vvect = interp3(XDFTvect,YDFTvect,ZDFTvect,VDFTsuper,Xvect,Yvect,Zvect,'linear');
    Vvect = InterpolatedFn(Xvect,Yvect,Zvect);%<- Beware! this step takes absolutely forever
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
  else
    %use new method

    N1 = linspace(0,12,Nxy);
    N2 = linspace(0,12,Nxy);
    N1_a = [N1 13];
    N2_a = [N2 13];

    %N1 = linspace(0,const.c,Nxy);
    %N2 = linspace(0,const.c,Nxy);
    %N1_a = [N1 const.c*(Nxy+1/Nxy)];
    %N2_a = [N2 const.c*(Nxy+1/Nxy)];
    [x1, y1, z1] =  ndgrid(N1_a, N1_a, Z);
    Vinterp_a = interp3(VDFT_a, x1, y1, z1, 'nearest');
    Vinterp = Vinterp_a(1:Nxy, 1:Nxy, :);
    figure
    equipotential_plot('V', Vinterp, 'V0', 0, 'z', Z, 'X', N1, 'Y', N2)
    shading interp
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
Ndefect = round(Nsites * Theta);
disp("Actual number of sites = " + Ndefect)
Nensemble = (factorial(Nsites)) ...
  /(factorial(Nsites - Ndefect)*factorial(Ndefect));
disp("Total ensemble size = " + Nensemble)
Nensemble_limit = 6;
if(Nensemble > Nensemble_limit)
  disp("Truncating ensemble to just " + Nensemble_limit)
  Nensemble = Nensemble_limit;
end
potStructArray = struct([]);
boolgrid_ensemble = zeros(Nsuper,Nsuper,Nensemble,'logical');
if(Ndefect == 0)
  potStructArray(1).V = Vsuper;
  potStructArray(1).a1=Nsuper*a1; potStructArray.a2=Nsuper*a2;
  potStructArray(1).zmin=Z(1);
  potStructArray(1).zmax=Z(end);
  potStructArray(1).zPoints=length(Z);

    %% Plot the potential
  plotPot = true;
  if(plotPot)
    Vplotted = Vsuper;
    %% Plot the potential. Disabled for now, as if the grid res is too high it complains
    %nPlot = 2/3;mPlot = 1/2;
    comparePots = true;
    nPlotDef = 0;mPlotDef = 0;
    aboveCol = [0.3 0. 1];

    nPlotHol = 2/3;mPlotHol = 1/3;
    holCol = [0.0 0.6 0.2];

    nPlotMo = 1/3;mPlotMo = 2/3;
    moCol = [1 0.2 0];

    if(comparePots)
      ComparePotentials(Vplotted,Vinterp,'Analytical potential','DFT interpolated',a1,a2,mPlotDef,nPlotDef,Z,Z,0,aboveCol)
      ComparePotentials(Vplotted,Vinterp,'Analytical potential','DFT interpolated',a1,a2,mPlotHol,nPlotHol,Z,Z,0,holCol)
      ComparePotentials(Vplotted,Vinterp,'Analytical potential','DFT interpolated',a1,a2,mPlotMo,nPlotMo,Z,Z,0,moCol)
    end
    % Plot of a slice of the potential in the nth row, that is for constant x
      row = floor(Nxy/2);
    figure
    contourf(Z,  linspace(0, const.c*Nsuper, Nxy*Nsuper), ...%!!!
        reshape(Vplotted(row,:,:), [Nxy*Nsuper,Nz]), linspace(-30,100,24))
    
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
      equipotential_plot('V', Vplotted, 'V0', Vsoup, 'z', Z, 'X', Xsuper, 'Y', Ysuper)
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
      
      %clf
    end
    %% Plot the potential
    fontsize(gcf,scale=1)
    zSample = 3.5;
    zRow = floor((zSample - zMin)/(zMax-zMin) * Nz);
    figure
    contourf(Xsuper,Ysuper,Vplotted(:,:,zRow),10)
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
    
  end
else
  %for the ensemble of potentials, how do we guarentee that they're
  %different???
  %we've got ms_ensemble and ns_ensemble, but how do we check wehther or
  %not they contain the same coordinates? I guec  ss what we do is take a
  %single set of the coordinates that we've generated, and then check them
  %against every single coordinate we have so far. Either that or create
  %some kind of weird information-theoretic checksum style thing to
  %automatically determine if two arrays of coordinate pairs are the same

  for Ne = 1:Nensemble
    solved = false;
    while(~solved)
      ms = squeeze(ones(1,Ndefect,1,1,'int8'))*69;
      ns = squeeze(ones(1,Ndefect,1,1,'int8'))*69;
      ms_available = 0:Nsuper-1;
      ns_available = 0:Nsuper-1;
      %Another question is - How do we decide which sites to add defects on to?
      %Do we directly implement avoidance of nearest-neighbors?
      avoidNearestNeighbors = false;
      if(avoidNearestNeighbors)
        error("Nearest neighbour avoidance not yet implemented!")
      end

      boolgrid = zeros(Nsuper,Nsuper,'logical');
      for d = 1:Ndefect
        disp("d = " + d)
        foundValidSpot = false;
        while(~foundValidSpot)
          m = randi(Nsuper)-1;
          n = randi(Nsuper)-1;
          if(boolgrid(m+1,n+1) == false)
            foundValidSpot = true;
          end
        end
        disp("m, n = ")
        disp([m n])
        ms(d) = m;
        ns(d) = n;
        %ms_available = ms_available(ms_available~=m);
        %ns_available = ns_available(ns_available~=n);
      end
      disp("Trial defects for ensemble " + Ne + " = ")
      disp([ms;ns])
      testgrid = zeros(Nsuper,Nsuper,'logical');
      for m = ms
        for n = ns
          testgrid(m+1,n+1) = true;
        end
      end
      solved = true;
      for ne = 1:Ne
        samplegrid = boolgrid_ensemble(:,:,Ne);
        if(samplegrid == testgrid)
          disp("This defect arrangement has already been stored!")
          solved = false;
        end
      end
    end
    disp("Solved!")
    boolgrid_ensemble(:,:,Ne) = testgrid;
  end
end
  
for Ne = 1:Nensemble
  disp("Ensemble number " + Ne)
  %Vsuper = AddSulphurDefect(false,Vsuper,1,1,a1,a2,Nsuper,Xsuper,Ysuper,Z);

  Vout = Vsuper;
  for m = 1:Nsuper
    for n = 1:Nsuper
      if(boolgrid_ensemble(m,n,Ne) == true)
        if((m == 0) || (n == 0) || (m == Nsuper-1) || (n == Nsuper-1))
          Vout = AddSulphurDefect(true,Vout,m,n,a1,a2,Nsuper,Xsuper,Ysuper,Z);
        else
          Vout = AddSulphurDefect(false,Vout,m,n,a1,a2,Nsuper,Xsuper,Ysuper,Z);
        end
      end
    end
  end
  potStructArray(Ne).V=Vout;
  potStructArray(Ne).a1=Nsuper*a1; potStructArray(Ne).a2=Nsuper*a2;
  potStructArray(Ne).zmin=Z(1);
  potStructArray(Ne).zmax=Z(end);
  potStructArray(Ne).zPoints=length(Z);
  %% Plot the potential
  plotPot = true;
  if(plotPot)
    Vplotted = Vout;
    %% Plot the potential. Disabled for now, as if the grid res is too high it complains
    %nPlot = 2/3;mPlot = 1/2;
    comparePots = false;
    nPlotDef = 0;mPlotDef = 0;
    aboveCol = [0.3 0. 1];

    nPlotHol = 2/3;mPlotHol = 1/3;
    holCol = [0.0 0.6 0.2];

    nPlotMo = 1/3;mPlotMo = 2/3;
    moCol = [1 0.2 0];

    if(comparePots)
      ComparePotentials(Vplotted,Vinterp,'Analytical potential','DFT interpolated',a1,a2,mPlotDef,nPlotDef,Z,Z,0,aboveCol)
      ComparePotentials(Vplotted,Vinterp,'Analytical potential','DFT interpolated',a1,a2,mPlotHol,nPlotHol,Z,Z,0,holCol)
      ComparePotentials(Vplotted,Vinterp,'Analytical potential','DFT interpolated',a1,a2,mPlotMo,nPlotMo,Z,Z,0,moCol)
    end
    % Plot of a slice of the potential in the nth row, that is for constant x
      row = floor(Nxy/2);
    figure
    contourf(Z,  linspace(0, const.c*Nsuper, Nxy*Nsuper), ...%!!!
        reshape(Vplotted(row,:,:), [Nxy*Nsuper,Nz]), linspace(-30,100,24))
    
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
      equipotential_plot('V', Vplotted, 'V0', Vsoup, 'z', Z, 'X', Xsuper, 'Y', Ysuper)
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
      
      %clf
    end
    %% Plot the potential
    fontsize(gcf,scale=1)
    zSample = 3;
    zRow = floor((zSample - zMin)/(zMax-zMin) * Nz);
    figure
    contourf(Xsuper,Ysuper,Vplotted(:,:,zRow),10)
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
    
  end
  %===
  
end
%===

%% for fitting
%Zdefect = dft.zAxis;
%Vdefect = dft.aboveDefect;
SpaghettiBolognaise = [a1(1) a2(1);a1(2) a2(2)]/(Nxy*Nsuper);
zFitMin = 1.5;
k = int8(interp1(Z,1:numel(Z),zFitMin));
if(k == 0)
  k = 1;
end
m = mPlotDef; n = nPlotDef;
centre = m*a1+n*a2;
result = SpaghettiBolognaise\(centre');
i = int8(result(1))+1;
j = int8(result(2))+1;
V1piece = squeeze(Vsuper(i,j,k:end));
Zpiece = Z(k:end);

weights = Zpiece;
for k = 1:numel(Zpiece)
  weights(k) = 1/(exp((3-Zpiece(k))*7)+1);
end
%figure
%plot(Zpiece,weights)
%% Get min and max bounds of the potentials
DFTmin = min(VDFTsuper,[],"all");
DFTmax = max(VDFTsuper,[],"all");
AnalyticMin = min(Vsuper,[],"all");
AnalyticMax = max(Vsuper,[],"all");

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
thetastr = [char(num2str(Theta))];
nens = [char(num2str(Nensemble))];
S = fileread('latticeVects.info_for_vivian_python_nice_plotting_hexagon_script');
realStr = ['Real space vectors:',newline,'a1 = ',a1str, newline, 'a2 = ',a2str,newline,'Nsuper = ',nsupstr];
recpStr = ['Reciprocal vectors:',newline,'b1 = ',b1str, newline, 'b2 = ', b2str];
defectStr = ['Defect data:',newline,'Theta = ', thetastr,newline,'Ensenble size = ', nens];
S = [realStr,newline,recpStr,newline,defectStr,S];
FID = fopen('latticeVects.info_for_vivian_python_nice_plotting_hexagon_script', 'w');
if FID == -1, error('Cannot open file %s', FileName); end
fwrite(FID, S, 'char');
fclose(FID);

%% We supply the lattice to the mulitscat script so it can do its thing
doingMSshit = true;
if(doingMSshit)
    %potStructArray.V = Vsuper;
    Multiscat.PreparePotentialFiles(potStructArray);
    
    Multiscat.prepareFourierLabels(Vsuper);
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
  function [V] = VSulph(z)
    D = 19.9886;
    a = 0.8122;
    alpha = 1.4477;
    b = 0.1958;
    beta = 0.2029;
    z0 = 3.3719;
    z1 = 1.7316;
    V = D*(exp(2*alpha*(z0-z))-2*a*exp(alpha*(z0-z))-2*b*exp(2*beta*(z1-z)));
  end

  function [V] = VHollow(z)
    D = 30.9674;
    a = 0.4641;
    alpha = 1.1029;
    b = 0.1993;
    beta = 0.6477;
    z0 = 3.1411;
    z1 = 3.8323;
    V = D*(exp(2*alpha*(z0-z))-2*a*exp(alpha*(z0-z))-2*b*exp(2*beta*(z1-z)));
  end

  function [V] = VMolyb(z)
    D = 30.8733;
    a = 0.9242;
    alpha = 1.1656;
    b = 0.0841;
    beta = 4.2710;
    z0 = 3.0991;
    z1 = 2.1265;
    V = D*(exp(2*alpha*(z0-z))-2*a*exp(alpha*(z0-z))-2*b*exp(2*beta*(z1-z)));
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
fittingDFT = true;
if(fittingDFT)
        %+ V1func(Z) * Qfunc(X,Y)...
    VmatrixElement = VSulph(Z) ... %blue, sulphur
       * Qhexfunc(X,Y) ...
       + VHollow(Z) ... %red, molybdenum
      * Qhexfunc(X,Y - (const.c/sqrt(3))) ...
      + VMolyb(Z) ...%green, hollow site
      * Qhexfunc(X-const.c/2,Y-(const.c*1/(2*sqrt(3))));
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

function [Vout] = AddSulphurDefect(doWeRepeat,Vin,min,nin,a1,a2,Nsuper,Xsuper,Ysuper,Z)
%Adds a defect at sulphur site (m,n)
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
  centre0 = double(min)*a1+double(nin)*a2;
  
  if doWeRepeat
    disp("Repeating!")
    for m = -1:1
      for n = -1:1
        centresX(m+2,n+2) = centre0(1)+m*a1(1)*Nsuper+n*a2(1)*Nsuper;
        centresY(m+2,n+2) = centre0(2)+m*a1(2)*Nsuper+n*a2(2)*Nsuper;
        for k = 1:Nz
          for i = 1:NxySuper
            for j = 1:NxySuper
              x = Xsuper(i,j);
              y = Ysuper(i,j);
              centre = [centresX(m+2,n+2) centresY(m+2,n+2)];
              Vout(i,j,k) = Vout(i,j,k)+val(x,y,Z,k,centre);
              %disp("x, y, z = " + x + ", " + y + ", " + Z(k) +...
              %     ", Value = " + val(x,y,Z,k,centre));
            end
          end
        end
      %disp("m, n = " + m + ", " + n)
      %disp("Centre:")
      %disp(centre)
      %equipotential_plot('V', Vout, 'V0', 0, 'z', Z, 'X', Xsuper, 'Y', Ysuper)
      %hold on
      %plot3(centresX(m+2,n+2),centresY(m+2,n+2),5,'.');
      %hold off
      %xlim([-3.5 2]*Nsuper);
      %ylim([-0.5 3]*Nsuper);
      %daspect([1 1 1])
      %shading interp
      end
    end
  else
    disp("Not repeating!")
    for k = 1:Nz
      for i = 1:NxySuper
        for j = 1:NxySuper
          x = Xsuper(i,j);
          y = Ysuper(i,j);
          Vout(i,j,k) = Vout(i,j,k)+val(x,y,Z,k,centre0);
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
  k1 = int8(interp1(Z1,1:numel(Z1),zMin))+1;
  centre = m*a1+n*a2;
  result = SpaghettiBolognaise1\(centre');
  i1 = int8(result(1))+1;
  j1 = int8(result(2))+1;
  V1piece = squeeze(V1(i1,j1,:));
  
  NxyNsuper2 = size(V2,1);
  SpaghettiBolognaise2 = [a1(1) a2(1);a1(2) a2(2)]/NxyNsuper2;
  k2 = int8(interp1(Z2,1:numel(Z2),zMin))+1;
  centre = m*a1+n*a2;
  result = SpaghettiBolognaise2\(centre');
  i2 = int8(result(1))+1;
  j2 = int8(result(2))+1;
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

function [checksum] = vivCheckSum(intArray1, intArray2)
  %Calculates checksum for two, non-ordered, but between them paired,
  %arrays of integers
  if(size(intArray1) ~= size(intArray2))
    error("Sizes of arrays do not match!")
  end
  length = size(intArray1,1);
  for i = 1:length
    m = intArray1(i);
    n = intArray2(i);
  end
end
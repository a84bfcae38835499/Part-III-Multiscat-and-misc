clear; close all; clc;
rng default;
rng("shuffle");
%Number of grid points, number of Z points, and number of lattices
%contained in the overall superlattice (or rather the square root of that)
Nxy = 32; Nz = 100; Nsuper = 2;
%Theta = 0.9;
Theta = (1/(Nsuper*Nsuper));
disp('Theta = ' + Theta)
Nensemble_limit = 1;
usingDisplacementDefects = true;
displacementMode = 1; % 0 = Gaussians
                      % 1 = Hemisphere
  defectH = 1.;
  defectR = 1.;
  minDist = defectR;
zMax = 6; zMin = 0;%units Å
fileprefix = "2x2sphere-test"
onlyWriteLatticeFile = false;
plotPot = true;
onlyPrepConf = false;

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

%% alskdjhshfjgd

nPlotDef = 0;mPlotDef = 0;
aboveCol = [0.3 0. 1];
nPlotHol = 2/3;mPlotHol = 1/3;
holCol = [0.0 0.6 0.2];
nPlotMo = 1/3;mPlotMo = 2/3;
moCol = [1 0.2 0];
nMidHol = 2/6;mMidHol = 1/6;
holMCol = [0.0 0.7 0.9];
nMidMo = 1/6;mMidMo = 2/6;
moMCol = [0.8 0.5 0];

%% Defect density calculations
cellArea = const.c^2 * sqrt(1-(dot(a1,a2)/(const.c^2))^2);
disp("Unit cell area = " + cellArea + "Å^2")
cellAreaS = cellArea * (Nsuper^2);
disp("Supercell area = " + cellAreaS + "Å^2")

disp("Target number of sites = " + (Nsuper*Nsuper * Theta))
Ndefect = int64(round(Nsuper*Nsuper * Theta));
disp("Actual number of sites = " + Ndefect)
defectDensity = double(Ndefect)/cellAreaS;
disp("defectDensity = " + defectDensity + "/Å^2")
defectDensity = defectDensity * ((1e10/1e2)^2);
disp("              = " + num2str(defectDensity,'%e') + "/cm^2")
disp("Theta = " + Theta)

probV = double(Ndefect)/double(Nsuper*Nsuper);
probS = 1 - probV;
if(probV > 0. && probS > 0.)
  inputEntropy = - double(Ndefect) * probV*log(probV) - ...
    (Nsuper^2-double(Ndefect)) * probS*log(probS);
  inputEntropy = inputEntropy/(-Nsuper^2*0.5*log(0.5));
else
  inputEntropy = 0;
end
disp("Input entropy = " + num2str(inputEntropy))
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

XDFTsuper = XDFTsuper - const.c/(sqrt(3)); %makes the 0,0 point be a sulphur
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

%% Now create pristine potentials
V = zeros(Nxy,Nxy,Nz);
Vinterp = zeros(Nxy,Nxy,Nz);
X = zeros(Nxy,Nxy);
Y = zeros(Nxy,Nxy);
Xsuper = zeros(Nxy*Nsuper,Nxy*Nsuper);
Ysuper = zeros(Nxy*Nsuper,Nxy*Nsuper);
Vsuper = zeros(Nsuper*Nxy,Nsuper*Nxy,Nz);
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

if(~usingDisplacementDefects)
  for k = 1:Nz
        V(:,:,k) = Vfunc_MoS2(X,Y,Z(k));
  end
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
end
%We strictly ought to be careful with boundary conditions cos MS doesn't
%actually check them lol
%===
%% Now interpolate the DFT data into a useful basis
interpolateDFTdata = false;
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
%Assume that we always have a defect at the (0,0) position, to fix
%translational invariance leadings to degeneracy
%if(Ndefect ~= 0 && (Nsuper*Nsuper)-1 - Ndefect > 0)
%  Nensemble = (factorial(Nsites-1)) ...
%    /(factorial(Nsites - Ndefect)*factorial(Ndefect));
  Nensemble = 10;  %gansta maths
%else
%  Nensemble = 8;
%end
disp("Total ensemble size = " + Nensemble)
if(Nensemble > Nensemble_limit)
  disp("Truncating ensemble to just " + Nensemble_limit)
  Nensemble = Nensemble_limit;
end



potStructArray = struct([]);
boolgrid_ensemble = zeros(Nsuper,Nsuper,Nensemble,'logical');
if(Ndefect == 0 || usingDisplacementDefects)
  %% 0 defects
  if(~usingDisplacementDefects)
    disp("No defects!!!")
    potStructArray(1).V = Vsuper;
    potStructArray(1).a1=Nsuper*a1; potStructArray.a2=Nsuper*a2;
    potStructArray(1).zmin=Z(1);
    potStructArray(1).zmax=Z(end);
    potStructArray(1).zPoints=length(Z);
    potStructArray(1).fileprefix=fileprefix;
    potStructArray(1).Nxy=Nxy;
    potStructArray(1).Nsuper=Nsuper;
    potStructArray(1).Ndefect=Ndefect;
  else
    disp("Using displacement defects!!!!")
    for Ne = 1:Nensemble
      randomX = 1337*ones(1,Ndefect);
      randomY = 1337*ones(1,Ndefect);
      repeat = zeros(1,Ndefect,'logical');
      for d = 1:Ndefect
        solved = false;
        while(~solved)
          r1 = rand*Nsuper;
          r2 = rand*Nsuper;
          if(r1 < defectR*2 || r1-defectR*2 > Nsuper*const.c || r2 < defectR*2 || r2-defectR*2 > Nsuper*const.c)
            repeat(d) = true;
          end
          R1 = r1 * a1;
          R2 = r2 * a2; 
          x = R1(1) + R2(1);
          y = R1(2) + R2(2);
          solved = true;
          for b = 1:d-1
            if(repeat(d))
                disp("Repeating!")
                for m = -1:1
                  for n = -1:1
                    centre = [randomX(b)+m*a1(1)*Nsuper+n*a2(1)*Nsuper
                            randomY(b)+m*a1(2)*Nsuper+n*a2(2)*Nsuper];
                    xt = centre(1);
                    yt = centre(2);
                    dist = (xt-x)^2+(yt-y)^2;
                    dist = sqrt(dist);
                    disp("dist = " + dist)
                    if(dist < minDist)
                      solved = false;
                    end
                  end
                end
              else
              xt = randomX(b);
              yt = randomY(b);
              dist = (xt-x)^2+(yt-y)^2;
              dist = sqrt(dist);
              disp("dist = " + dist)
              if(dist < minDist)
                solved = false;
              end
            end
          end
        disp("solved = " + solved);
        end
        disp("Solved!");
        randomX(d) = x;
        randomY(d) = y;
      end
      addZ = zeros(Nxy*Nsuper,Nxy*Nsuper);
      for d = 1:Ndefect
        if(repeat(d))
          disp("Repeating!")
          for m = -1:1
            for n = -1:1
              centre = [randomX(d)+m*a1(1)*Nsuper+n*a2(1)*Nsuper
                      randomY(d)+m*a1(2)*Nsuper+n*a2(2)*Nsuper];
              if(displacementMode == 0)
                    addZ = addZ - defectH*Gaussian2D(Xsuper,Ysuper,centre,defectR*3);
              elseif(displacementMode == 1)
                  ikbT = 10;
                  mu = defectR*.75;
                  xt = centre(1);
                  yt = centre(2);
                  dist = (xt-Xsuper).^2+(yt-Ysuper).^2;
                  r = sqrt(dist);
                  yarr = max(0,1-dist/(defectR^2));
                  yarr = sqrt(yarr);

                  factor = (1./( 1+exp((r-mu)*ikbT) ));
                  factor = (factor.*( 1+exp((-mu)*ikbT) ));
                  yarr = factor.*(yarr)+ ... 
                    (1-factor).*(Gaussian2D(Xsuper,Ysuper,centre,defectR/2));
                  addZ = addZ - defectH*yarr;
              end
               
            end
          end
        else
          centre = [randomX(d)
                  randomY(d)];
          if(displacementMode == 0)
                addZ = addZ - defectH*Gaussian2D(Xsuper,Ysuper,centre,defectR/2);
          elseif(displacementMode == 1)
                  ikbT = 10;
                  mu = defectR*.75;
                  xt = centre(1);
                  yt = centre(2);
                  dist = (xt-Xsuper).^2+(yt-Ysuper).^2;
                  r = sqrt(dist);
                  yarr = max(0,1-dist/(defectR^2));
                  yarr = sqrt(yarr);

                  factor = (1./( 1+exp((r-mu)*ikbT) ));
                  factor = (factor.*( 1+exp((-mu)*ikbT) ));
                  yarr = factor.*(yarr)+ ... 
                    (1-factor).*(Gaussian2D(Xsuper,Ysuper,centre,defectR/2));
                  addZ = addZ - defectH*yarr;
          end
        end
      end
      for i = 1:Nxy*Nsuper
        for j = 1:Nxy*Nsuper
          for k = 1:Nz
            Vsuper(i,j,k) = Vfunc_MoS2(Xsuper(i,j),Ysuper(i,j),Z(k) + addZ(i,j));
          end
        end
      end
      potStructArray(Ne).V = Vsuper;
      potStructArray(Ne).a1=Nsuper*a1; potStructArray(Ne).a2=Nsuper*a2;
      potStructArray(Ne).zmin=Z(1);
      potStructArray(Ne).zmax=Z(end);
      potStructArray(Ne).zPoints=length(Z);
      potStructArray(Ne).fileprefix=fileprefix;
      potStructArray(Ne).Nxy=Nxy;
      potStructArray(Ne).Nsuper=Nsuper;
      potStructArray(Ne).Ndefect=Ndefect;
  if(plotPot)
    Vplotted = Vsuper;
    %nPlot = 2/3;mPlot = 1/2;
    comparePots = false;
    if(comparePots)
      [xS, yS] = ComparePotentials(Vplotted,dft.aboveSd,'Analytical potential','DFT - Vacancy',a1,a2,mPlotDef,nPlotDef,Z,dft.zAxis,0,aboveCol,Nxy);
      [xH, yH] = ComparePotentials(Vplotted,dft.aboveHollowd,'Analytical potential','DFT - Hollow site',a1,a2,mPlotHol,nPlotHol,Z,dft.zAxis,0,holCol,Nxy);
      [xM, yM] = ComparePotentials(Vplotted,dft.aboveMod,'Analytical potential','DFT - Molybdenum',a1,a2,mPlotMo,nPlotMo,Z,dft.zAxis,0,moCol,Nxy);
      [xHm, yHm] = ComparePotentials(Vplotted,dft.midHo,'Analytical potential','DFT - Mid Hollow site',a1,a2,mMidHol,nMidHol,Z,dft.zAxisHiRes,0,holMCol,Nxy);
      [xMm, yMm] = ComparePotentials(Vplotted,dft.midMo,'Analytical potential','DFT - Mid Molybdenum',a1,a2,mMidMo,nMidMo,Z,dft.zAxisHiRes,0,moMCol,Nxy);
    end
    % Plot of a slice of the potential in the nth row, that is for constant x
      row = floor(Nxy/2);
    figure
    contourf(Z,  linspace(0, const.c*Nsuper, Nxy*Nsuper), ...%!!!
        reshape(Vplotted(row,:,:), [Nxy*Nsuper,Nz]), linspace(-30,100,24))
    
        fontsize(gcf,scale=1)
    xlabel('z/Å')
    ylabel('y/Å') %is this x or y? I think y but idrk;
    colorbar
    xlim([1.5,6])
    title('Potential in z, used in simulation')
    hbar = colorbar;
    ylabel(hbar,'Energy / meV');
    figure
    fileindx = 1;
    for i = 1e-5
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
      savestr = "Figures/Frames/frame_" +num2str(fileindx,'%06d')+ ".jpg";
      fileindx = fileindx + 1;
      saveas(gcf,savestr,'jpg')
    end
    fontsize(gcf,scale=1)
    zSample = 2.5;
    zRow = floor((zSample - zMin)/(zMax-zMin) * Nz);
    figure
    contourf(Xsuper,Ysuper,Vplotted(:,:,zRow),16)
    daspect([1 1 1])
    xlabel('x/Å')
    ylabel('y/Å')
    title('Potentials at z = ' + string(zSample) + ' Å');
    colormap(parula(16))
    hbar = colorbar('southoutside');
    xlabel(hbar,'Energy / meV');
    %add indicators for where we're sampling the potential z
    fontsize(gcf,scale=1)

    plotPoints = false;
    if(plotPoints)
      hold on
      %xPlot = mPlotDef*a1(1)+nPlotDef*a2(1);
      %yPlot = mPlotDef*a1(2)+nPlotDef*a2(2);
      %plot(xPlot,yPlot,'*',MarkerSize=24,Color=aboveCol);
      %plot(xPlot,yPlot,'.',MarkerSize=24,Color=aboveCol);
      plot(xS,yS,'*',MarkerSize=24,Color=aboveCol);
      plot(xS,yS,'.',MarkerSize=24,Color=aboveCol);
    
      %xPlot = mPlotHol*a1(1)+nPlotHol*a2(1);
      %yPlot = mPlotHol*a1(2)+nPlotHol*a2(2);
      %plot(xPlot,yPlot,'*',MarkerSize=24,Color=holCol);
      %plot(xPlot,yPlot,'.',MarkerSize=24,Color=holCol);
      plot(xH,yH,'+',MarkerSize=24,Color=holCol);
      plot(xH,yH,'.',MarkerSize=24,Color=holCol);
    
      %xPlot = mPlotMo*a1(1)+nPlotMo*a2(1);
      %yPlot = mPlotMo*a1(2)+nPlotMo*a2(2);
      %plot(xPlot,yPlot,'*',MarkerSize=24,Color=moCol);
      %plot(xPlot,yPlot,'.',MarkerSize=24,Color=moCol);
      plot(xM,yM,'x',MarkerSize=24,Color=moCol);
      plot(xM,yM,'.',MarkerSize=24,Color=moCol);

      %xPlot = mMidHol*a1(1)+nMidHol*a2(1);
      %yPlot = mMidHol*a1(2)+nMidHol*a2(2);
      %plot(xPlot,yPlot,'*',MarkerSize=24,Color=holMCol);
      %plot(xPlot,yPlot,'.',MarkerSize=24,Color=holMCol);
      plot(xHm,yHm,'d',MarkerSize=24,Color=holMCol);
      plot(xHm,yHm,'.',MarkerSize=24,Color=holMCol);

      %xPlot = mMidMo*a1(1)+nMidMo*a2(1);
      %yPlot = mMidMo*a1(2)+nMidMo*a2(2);
      %plot(xPlot,yPlot,'*',MarkerSize=24,Color=moMCol);
      %plot(xPlot,yPlot,'.',MarkerSize=24,Color=moMCol);
      plot(xMm,yMm,'p',MarkerSize=24,Color=moMCol);
      plot(xMm,yMm,'.',MarkerSize=24,Color=moMCol);

      defCol = [0.3 0.7 1];
      for ne = 1:Nensemble
        for m = 0:int64(Nsuper-1)
          for n = 0:int64(Nsuper-1)
            if(boolgrid_ensemble(m+1,n+1,Ne))
              xPlot = double(m)*a1(1)+double(n)*a2(1);
              yPlot = double(m)*a1(2)+double(n)*a2(2);
              plot(xPlot,yPlot,'p',MarkerSize=24, MarkerEdgeColor=[0 0 0],MarkerFaceColor=defCol);
            end
          end
        end
      end
      hold off

    end
  savestr = "Figures/" + fileprefix + "_" + string(Ne) + ".jpg";
  saveas(gcf,savestr,'jpg')
  end
    end
  end
else
  %% One or more defects
  disp("One or more defects!!!")
  %for the ensemble of potentials, how do we guarentee that they're
  %different???
  %we've got ms_ensemble and ns_ensemble, but how do we check wehther or
  %not they contain the same coordinates? I guec  ss what we do is take a
  %single set of the coordinates that we've generated, and then check them
  %against every single coordinate we have so far. Either that or create
  %some kind of weird information-theoretic checksum style thing to
  %automatically determine if two arrays of coordinate pairs are the same
  successful = 0;
  for Ne = 1:Nensemble
    solved = false;
    pizza = 0;
    while(~solved)
      pizza = pizza + 1;
      if(pizza > 500)
        disp("pizza = " + pizza)
        Nensemble = successful;
        error("New Nensemble = " + successful)
        solved = true;
      else
      ms = squeeze(ones(Ndefect,1,1,1,'int64'))*69;
      ns = squeeze(ones(Ndefect,1,1,1,'int64'))*69;
      ms_available = 0:Nsuper-1;
      ns_available = 0:Nsuper-1;
      %Another question is - How do we decide which sites to add defects on to?
      %Do we directly implement avoidance of nearest-neighbors?
      
      boolgrid = zeros(Nsuper,Nsuper,'logical');
      
      boolgrid(1,1) = true;
      ms(1) = 0;
      ns(1) = 0;
    
          for d = 2:Ndefect
            %disp("d = " + d)
            foundValidSpot = false;
            while(~foundValidSpot)
              m = randi(Nsuper)-1;
              n = randi(Nsuper)-1;
              avoidNearestNeighbors = true;
              if(avoidNearestNeighbors)
                if(boolgrid(m+1,n+1) == false ...
                    && boolgrid(mod(m+1,Nsuper)+1,n+1) == false ...
                    && boolgrid(mod(m-1,Nsuper)+1,n+1) == false ...
                    && boolgrid(mod(m,Nsuper)+1,n+1) == false ...
                    && boolgrid(m+1,mod(n+1,Nsuper)+1) == false ...
                    && boolgrid(m+1,mod(n-1,Nsuper)+1) == false ...
                    && boolgrid(mod(m+1,Nsuper)+1,mod(n+1,Nsuper)+1) == false ...
                    && boolgrid(mod(m-1,Nsuper)+1,mod(n-1,Nsuper)+1) == false)
                  foundValidSpot = true;
                  boolgrid(m+1,n+1) = true;
                end
              else
                if(boolgrid(m+1,n+1) == false)
                  foundValidSpot = true;
                  boolgrid(m+1,n+1) = true;
                end
              end
            end
            %disp("m, n = ")
            %disp([m n])
            ms(d) = m;
            ns(d) = n;
            %ms_available = ms_available(ms_available~=m);
            %ns_available = ns_available(ns_available~=n);
          end
          
          testgrid = zeros(Nsuper,Nsuper,'logical');
          for index = 1:length(ms)
            mt = ms(index)+1;
            nt = ns(index)+1;
            testgrid(mt,nt) = true;
          end
          %disp("Trial defect arrangement for ensemble " + Ne + " = ")
          %disp(testgrid)
          solved = true;
          for ne = 1:Ne
            samplegrid = boolgrid_ensemble(:,:,ne);
            %disp("Samplegrid = ")
            %disp(samplegrid)
            if( AreCyclicBoundaryMatriciesEqual(samplegrid,testgrid,Ndefect))
              disp("This defect arrangement has already been stored!")
              solved = false;
            end
          end
          disp("Solved = " + solved)
      end
      end
        disp("Solved!")
        successful = successful + 1;
        boolgrid_ensemble(:,:,Ne) = testgrid;
    end
for Ne = 1:Nensemble
  disp("Ensemble number " + Ne)
  %Vsuper = AddSulphurDefect(false,Vsuper,1,1,a1,a2,Nsuper,Xsuper,Ysuper,Z);

  Vout = Vsuper;
  for m = 0:int64(Nsuper-1)
    for n = 0:int64(Nsuper-1)
      if(boolgrid_ensemble(m+1,n+1,Ne) == true)
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
  potStructArray(Ne).Nxy=Nxy;
  potStructArray(Ne).Nsuper=Nsuper;
  potStructArray(Ne).Ndefect=Ndefect;
  potStructArray(Ne).fileprefix=fileprefix;
  if(plotPot)
    Vplotted = Vout;
    %nPlot = 2/3;mPlot = 1/2;
    comparePots = false;
    if(comparePots)
      [xS, yS] = ComparePotentials(Vplotted,dft.aboveSd,'Analytical potential','DFT - Vacancy',a1,a2,mPlotDef,nPlotDef,Z,dft.zAxis,0,aboveCol,Nxy);
      [xH, yH] = ComparePotentials(Vplotted,dft.aboveHollowd,'Analytical potential','DFT - Hollow site',a1,a2,mPlotHol,nPlotHol,Z,dft.zAxis,0,holCol,Nxy);
      [xM, yM] = ComparePotentials(Vplotted,dft.aboveMod,'Analytical potential','DFT - Molybdenum',a1,a2,mPlotMo,nPlotMo,Z,dft.zAxis,0,moCol,Nxy);
      [xHm, yHm] = ComparePotentials(Vplotted,dft.midHo,'Analytical potential','DFT - Mid Hollow site',a1,a2,mMidHol,nMidHol,Z,dft.zAxisHiRes,0,holMCol,Nxy);
      [xMm, yMm] = ComparePotentials(Vplotted,dft.midMo,'Analytical potential','DFT - Mid Molybdenum',a1,a2,mMidMo,nMidMo,Z,dft.zAxisHiRes,0,moMCol,Nxy);
    end
    % Plot of a slice of the potential in the nth row, that is for constant x
      row = floor(Nxy/2);
    figure
    contourf(Z,  linspace(0, const.c*Nsuper, Nxy*Nsuper), ...%!!!
        reshape(Vplotted(row,:,:), [Nxy*Nsuper,Nz]), linspace(-30,100,24))
    
        fontsize(gcf,scale=1)
    xlabel('z/Å')
    ylabel('y/Å') %is this x or y? I think y but idrk;
    colorbar
    xlim([1.5,6])
    title('Potential in z, used in simulation')
    hbar = colorbar;
    ylabel(hbar,'Energy / meV');
    figure
    fileindx = 1;
    for i = 1e-5
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
      savestr = "Figures/Frames/frame_" +num2str(fileindx,'%06d')+ ".jpg";
      fileindx = fileindx + 1;
      saveas(gcf,savestr,'jpg')
    end
    fontsize(gcf,scale=1)
    zSample = 2.5;
    zRow = floor((zSample - zMin)/(zMax-zMin) * Nz);
    figure
    contourf(Xsuper,Ysuper,Vplotted(:,:,zRow),16)
    daspect([1 1 1])
    xlabel('x/Å')
    ylabel('y/Å')
    title('Potentials at z = ' + string(zSample) + ' Å');
    colormap(parula(16))
    hbar = colorbar('southoutside');
    xlabel(hbar,'Energy / meV');
    %add indicators for where we're sampling the potential z
    fontsize(gcf,scale=1)

    plotPoints = false;
    if(plotPoints)
      hold on
      %xPlot = mPlotDef*a1(1)+nPlotDef*a2(1);
      %yPlot = mPlotDef*a1(2)+nPlotDef*a2(2);
      %plot(xPlot,yPlot,'*',MarkerSize=24,Color=aboveCol);
      %plot(xPlot,yPlot,'.',MarkerSize=24,Color=aboveCol);
      plot(xS,yS,'*',MarkerSize=24,Color=aboveCol);
      plot(xS,yS,'.',MarkerSize=24,Color=aboveCol);
    
      %xPlot = mPlotHol*a1(1)+nPlotHol*a2(1);
      %yPlot = mPlotHol*a1(2)+nPlotHol*a2(2);
      %plot(xPlot,yPlot,'*',MarkerSize=24,Color=holCol);
      %plot(xPlot,yPlot,'.',MarkerSize=24,Color=holCol);
      plot(xH,yH,'+',MarkerSize=24,Color=holCol);
      plot(xH,yH,'.',MarkerSize=24,Color=holCol);
    
      %xPlot = mPlotMo*a1(1)+nPlotMo*a2(1);
      %yPlot = mPlotMo*a1(2)+nPlotMo*a2(2);
      %plot(xPlot,yPlot,'*',MarkerSize=24,Color=moCol);
      %plot(xPlot,yPlot,'.',MarkerSize=24,Color=moCol);
      plot(xM,yM,'x',MarkerSize=24,Color=moCol);
      plot(xM,yM,'.',MarkerSize=24,Color=moCol);

      %xPlot = mMidHol*a1(1)+nMidHol*a2(1);
      %yPlot = mMidHol*a1(2)+nMidHol*a2(2);
      %plot(xPlot,yPlot,'*',MarkerSize=24,Color=holMCol);
      %plot(xPlot,yPlot,'.',MarkerSize=24,Color=holMCol);
      plot(xHm,yHm,'d',MarkerSize=24,Color=holMCol);
      plot(xHm,yHm,'.',MarkerSize=24,Color=holMCol);

      %xPlot = mMidMo*a1(1)+nMidMo*a2(1);
      %yPlot = mMidMo*a1(2)+nMidMo*a2(2);
      %plot(xPlot,yPlot,'*',MarkerSize=24,Color=moMCol);
      %plot(xPlot,yPlot,'.',MarkerSize=24,Color=moMCol);
      plot(xMm,yMm,'p',MarkerSize=24,Color=moMCol);
      plot(xMm,yMm,'.',MarkerSize=24,Color=moMCol);

      defCol = [0.3 0.7 1];
      for ne = 1:Nensemble
        for m = 0:int64(Nsuper-1)
          for n = 0:int64(Nsuper-1)
            if(boolgrid_ensemble(m+1,n+1,Ne))
              xPlot = double(m)*a1(1)+double(n)*a2(1);
              yPlot = double(m)*a1(2)+double(n)*a2(2);
              plot(xPlot,yPlot,'p',MarkerSize=24, MarkerEdgeColor=[0 0 0],MarkerFaceColor=defCol);
            end
          end
        end
      end
      hold off

    end
  savestr = "Figures/" + fileprefix + "_" + string(Ne) + ".jpg";
  saveas(gcf,savestr,'jpg')
  end
  %===
  
end
%===

end
%% data for python hex plotter
WritePythonInfo(fileprefix,a1,a2,cellArea,b1,b2,Nsuper,Theta,Nensemble,inputEntropy,defectDensity,Ndefect);
if(onlyWriteLatticeFile)
    error("Done!")
end
  
%% for fitting
%Zdefect = dft.zAxis;
%Vdefect = dft.aboveDefect;

%nPlot = 2/3;mPlot = 1/2;
comparePots = true;
SpaghettiBolognaise = [a1(1) a2(1);a1(2) a2(2)]/(Nxy*Nsuper);
zFitMin = 1.5;
k = int64(interp1(Z,1:numel(Z),zFitMin));
if(k == 0)
  k = 1;
end
m = mPlotDef; n = nPlotDef;
centre = m*a1+n*a2;
result = SpaghettiBolognaise\(centre');
i = int64(result(1))+1;
j = int64(result(2))+1;
%V1piece = squeeze(Vsuper(i,j,k:end));
%Zpiece = Z(k:end);
Zpiece = dft.zAxis;
VDefectpiece = dft.aboveMod;
VPristinepiece = dft.aboveMop;
weights = Zpiece;
for k = 1:numel(Zpiece)
  weights(k) = 1/(exp((2-Zpiece(k))*6)+1);
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
%% Plot corrugation
%[newx,newy] = meshgrid(-const.c:0.01:const.c,-const.c:0.01:const.c);
%newz = Qhexfunc(newx,newy);
%r = sqrt(newx.^2+newy.^2)/const.c;
%newz = newz .* (1./(1+exp((1/sqrt(3)-r)*1000)));
%surf(newx,newy,newz)
%shading interp
%colormap(plasma)
%daspect([1 1 1])
%% We supply the lattice to the mulitscat script so it can do its thing
doingMSshit = true;
if(doingMSshit)
    %potStructArray.V = Vsuper;
    confStruct=Multiscat.createConfigStruct(potStructArray);
    Multiscat.prepareConfigFile(confStruct);
    if(onlyPrepConf) 
        error("Done!")
    end
    Multiscat.prepareFourierLabels(Vsuper,fileprefix);
    Multiscat.PreparePotentialFiles(potStructArray);
end
%===
%% Function definitions

function [b1,b2,b3] = Reciprocal(a1,a2,a3)
    factor = 2*pi/dot(cross(a1,a2),a3);
    b1 = factor*cross(a2,a3);
    b2 = factor*cross(a3,a1);
    b3 = factor*cross(a1,a2);
end

function [VmatrixElement] = Vfunc_LiF(x,y,z)
    function [V0] = V0func(z)
        V0 = const.D * exp(2*const.alpha*(const.z0-z))...
            - 2*const.D*exp(const.alpha*(const.z0-z));
    end
    function [V1] = V1func(z)
        V1 = -2*const.beta*const.D*exp(2*const.alpha*(const.z0-z));
    end
    function [Q] = Qfunc(x,y)
        Q = cos(2*pi*x/const.a) + cos(2*pi*y/const.a);
    end
    VmatrixElement = V0func(z) + V1func(z)...
        * Qfunc(x,y);
end

function [VmatrixElement] = Vfunc_MoS2(X,Y,Z)
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
    D = 24.9674;
    a = 0.4641;
    alpha = 1.1029;
    b = 0.1993;
    beta = 0.6477;
    z0 = 3.1411;
    z1 = 3.8323;
    V = D*(exp(2*alpha*(z0-z))-2*a*exp(alpha*(z0-z))-2*b*exp(2*beta*(z1-z)));
  end

  function [V] = VMolyb(z)
    D = 20.1000;
    a = 0.9996;
    alpha = 1.1500;
    b = 0.0026;
    beta = 1.2439;
    z0 = 3.2200;
    z1 = 4.1864;
    V = D*(exp(2*alpha*(z0-z))-2*a*exp(alpha*(z0-z))-2*b*exp(2*beta*(z1-z)));
  end
        %+ V1func(Z) * Qfunc(X,Y)...
    VmatrixElement = VSulph(Z) ... %blue, sulphur
       * Qhexfunc(X,Y) ...
       + VHollow(Z) ... %green, hollow site
      * Qhexfunc(X,Y - (const.c/sqrt(3))) ...
      + VMolyb(Z) ...%red, molybdenum
      * Qhexfunc(X-const.c/2,Y-(const.c*1/(2*sqrt(3))));
end

function [Vout] = AddSulphurDefect(doWeRepeat,Vin,m_in,n_in,a1,a2,Nsuper,Xsuper,Ysuper,Z)
%Adds a defect at sulphur site (m,n)
  Vout = Vin;
  NxySuper = size(Vout,1);
  Nz = size(Vout,3);
  centre0 = double(m_in)*a1+double(n_in)*a2;
  if doWeRepeat
    disp("Repeating!")
    for m = -1:1
      for n = -1:1
        for k = 1:Nz
          centre = [centre0(1)+m*a1(1)*Nsuper+n*a2(1)*Nsuper
                    centre0(2)+m*a1(2)*Nsuper+n*a2(2)*Nsuper];
          %r = (x-centre(1))^2+(y-centre(2))^2;
          %r = sqrt(r)/const.c;
          %disp("r = " + num2str(r))
          %disp("vmatrixelem = " + Vout(i,j,k))
          %disp("defect val  = " + val(x,y,Z(k),centre))
          %disp("sum         = " + num2str(Vout(i,j,k) + val(x,y,Z(k),centre)));
          Vout(:,:,k) = Vout(:,:,k)+val(Xsuper,Ysuper,Z(k),centre);
          %disp("------------- ")
          %disp("x, y, z = " + x + ", " + y + ", " + Z(k) +...
          %     ", Value = " + val(x,y,Z,k,centre));
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
        Vout(:,:,k) = Vout(:,:,k)+val(Xsuper,Ysuper,Z(k),centre0);
        %disp("x, y, z = " + x + ", " + y + ", " + Z(k) +...
        %    ", Value = " + val);
    end
  end

  function [v] = val(x,y,z,centre)
    
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
    
    r = (x-centre(1)).^2+(y-centre(2)).^2;
    r = sqrt(r)./const.c;
    %maxr = min(r,[],"all");
    %minr = min(r,[],"all");
    extent = 0.3;
    cutoff = 1;
    s = extent * const.c;
    v = 0;
    c = 0.0311;
    d = 32.8260;
    e = 16.3770;
    gamma	= 0.9209;
    lambda = 0.9204;
    z2 = 5.6012;
    z3 = 3.7072;
    %r = max(0.1,r);

    flatDefect = true;
    hexDefect = false;
    if(flatDefect)
      if(hexDefect)
        ikbT = 12.9;
        mu = 0.92;
      else
        ikbT = 15.9;
        mu = 0.49;

        %ikbT = 4;
        %mu = 0.5;
      end
      VmatrixElement = Vfunc_MoS2(x,y,z);
    else
      error("This is impossible to fit, don't use this")
      if(hexDefect)
        ikbT = 17.9;
        mu = 1.6;
      else
        ikbT = 15.9;
        mu = 0.6;
      end
      VmatrixElement = VSulph(z) * Qhexfunc(x,y);
    end
    
    %args = (y-centre(2))./(x-centre(1));
    %angle = arrayfun(@(arg) atan(arg),args);
    %angle(isnan(angle))=0;
    %cutoffR = 1/sqrt(3)*cos(pi/6)./(cos(angle-(2*pi*floor((6*angle+pi)/(2*pi)))/6));
    %cutoffR = mu .* cutoffR;

    factor = (1./( 1+exp((r-mu)*ikbT) ));
    factor = (factor.*( 1+exp((-mu)*ikbT) ));
    %factor = factor./( 1+exp((r-hardCut)*10000));
    %maxf = max(factor,[],"all");
    %minf = min(factor,[],"all");
    %if(maxf > 1. )
    %  error("Max more than 1!")
    %elseif(minf < 0.)
    %    error("Min less than 0!")
    %end
    if(isnan(factor))
      error("Nans found!")
    end

    v = (-VmatrixElement + d*(exp(2*gamma*(z2-z))-2*c*exp(gamma*(z2-z)) ...
      -2*e*exp(2*lambda*(z3-z)))).*factor;
    %disp("Maxv = " + max(v,[],"all"))
    %disp("Minv = " + min(v,[],"all"))
    %disp("1./( 1+exp((-0.45)*10) = ")
    %disp(1./( 1+exp(-0.45*10) ))
  end
end

function [x,y] = ComparePotentials(Vgridded,Vline,V1name,V2name,a1,a2,m,n,Z1,Z2,zMin,plotColor,Nxy)
  Nsuper = double(size(Vgridded,1))/Nxy;
  SpaghettiBolognaise = [a1(1) a2(1);a1(2) a2(2)]/Nxy;
  k1 = int64(interp1(Z1,1:numel(Z1),zMin))+1;
  centre = m*a1+n*a2;
  result = SpaghettiBolognaise\(centre');
  i = int64(result(1));
  j = int64(result(2));
  V1piece = squeeze(Vgridded(mod(i+1,Nxy*Nsuper),mod(j+1,Nxy*Nsuper),:));
  x = double(i) * a1(1) + double(j) * a2(1);
  y = double(i) * a1(2) + double(j) * a2(2);
  x = x / Nxy;
  y = y / Nxy;
  x_ = centre(1);
  y_ = centre(2);
  disp("x, y = " + x + ", " + y)
  disp("x_,y_= " + x_ + ", " + y_)
  k2 = int64(interp1(Z2,1:numel(Z2),zMin))+1;
  %V2piece = squeeze(V2(i2,j2,:));
  V2piece = Vline;
  figure
  p = plot(Z1(k1:end),V1piece(k1:end),'DisplayName',V1name);
  p.LineStyle = ":";
  p.Color = plotColor;
  p.Marker = ".";
  p.LineWidth=2;
  tit = 'Plot of potentials at x, y = ' + string(x) + ', ' + string(y);
  title(tit)
  xlabel('z/Å')
  ylabel('Energy/meV')
  hold on
  legend()
  %disp(size(Z2))
  %disp(size(V2piece))
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
function [Q] = Qhexfunc(X,Y)
  %X_n = X ./ (const.c);
  %Y_n = Y ./ (const.c*sqrt(3));
  %Q = ((cos(2*pi*(X_n-Y_n))+cos(4*pi*Y_n)+cos(2*pi*(X_n+Y_n))) + 3/2)/(4.5);
  Q = 1;
end
clear; close all; clc;
rng(3);

%a = 2.84Å. see const.m for more stuff
a1=[const.a,0];
a2=[0,const.a];
%a1=[const.a,0];
%a2=[const.a/2,const.a * (sqrt(3)/2)];
a3=[0,0,const.a];
[b1,b2,b3] = Reciprocal([a1,0],[a2,0],a3);
%Number of grid points, number of Z points, and number of lattices
%contained in the overall superlattice (or rather the square root of that)
Nxy = 64; Nz = 100; Nsuper = 1;
zMax = 6; zMin = -2;%units Å

V = zeros(Nxy,Nxy,Nz);
X = zeros(Nxy,Nxy);
Y = zeros(Nxy,Nxy);
Xsuper = zeros(Nxy*Nsuper,Nxy*Nsuper);
Ysuper = zeros(Nxy*Nsuper,Nxy*Nsuper);
noiseX = dsp.ColoredNoise(0,NumChannels=Nxy*Nsuper,SamplesPerFrame=Nxy*Nsuper);
noiseY = dsp.ColoredNoise(0,NumChannels=Nxy*Nsuper,SamplesPerFrame=Nxy*Nsuper);
%noiseField = wgn(Nxy*Nsuper,Nxy*Nsuper,1,1,69420);%can be both positive and negative
noiseField = noiseX() + transpose(noiseY());

noiseFieldSuper = [noiseField noiseField noiseField; noiseField noiseField noiseField; noiseField noiseField noiseField];


noiseFieldSuper = smoothdata2(noiseFieldSuper,SmoothingFactor=0.25);


noiseField = noiseFieldSuper(Nxy*Nsuper:2*Nxy*Nsuper,Nxy*Nsuper:2*Nxy*Nsuper);
noiseMax = max(max(abs(noiseField)));
noiseField = 2*noiseField/noiseMax;


imagesc(noiseField);
%maxNoise = max(max(noiseField))
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

if false
  for k = 1:Nz
      V(:,:,k) = Vfunc(X,Y,Z(k));
  end
else
  for k = 1:Nz
        V(:,:,k) = Vnoise(X,Y,Z(k),noiseField);
  end
end
%We strictly ought to be careful with boundary conditions cos MS doesn't
%actually check them lol
%===
%Now we duplicate the lattice to get a superlattice

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
%Now add imperfections to the lattice
if false
  for k = 1:size(V,3) %Should be the z direction
    dropoff = Dropoff(Z(k));
    for i = 1:Nxy*Nsuper
      for j = 1:Nxy*Nsuper
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
%We supply the lattice to the mulitscat script so it can do its thing

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
%We also prepare a .csv which contains an equipotential plot.
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
%All this stuff is just to write a comment line lmao
writematrix(equipotentialMat,'Equipotential.csv','Delimiter', ',')
S = fileread('Equipotential.csv');
S = ['#Energy surface at ', eqCharArr, newline, S];
FID = fopen('Equipotential.csv', 'w');
if FID == -1, error('Cannot open file %s', FileName); end
fwrite(FID, S, 'char');
fclose(FID);

%===

%% Plot the potential
% Plot of a slice of the potential in the nth row
row = 1;
figure
contourf(Z,  linspace(0, const.a*Nsuper, Nxy*Nsuper), ...
    reshape(Vsuper(row,:,:), [Nxy*Nsuper,Nz]), linspace(-20,100,24))
xlabel('z/Å')
ylabel('x/Å')
colorbar
xlim([-1,4])
title('Potential in z, used in simulation')
hbar = colorbar;
ylabel(hbar,'Energy / meV');

    fontsize(gcf,scale=2)
% Linearly interpolated equipotential plot
equipotential_plot('V', V, 'z', Z, 'X', X, 'Y', Y)


%Function definitions

function [b1,b2,b3] = Reciprocal(a1,a2,a3)
    factor = 2*pi/dot(cross(a1,a2),a3);
    b1 = factor*cross(a2,a3);
    b2 = factor*cross(a3,a1);
    b3 = factor*cross(a1,a2);
end

function [VmatrixElement] = Vfunc(x,y,z)
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

function[Vmat] = Vnoise(X,Y,z,noiseArr)
      function [V0] = V0func(z)
        V0 = const.D * exp(2*const.alpha*(const.z0-z))...
            - 2*const.D*exp(const.alpha*(const.z0-z));
    end
    function [V1] = V1func(z)
        V1 = -2*const.beta*const.D*exp(2*const.alpha*(const.z0-z));
     
    end
  Vmat = zeros(size(X,1),size(Y,1));
  for i = 1:size(X,1)
    for j = 1:size(Y,1)
      Vmat(i,j) = V0func(z) + V1func(z) * noiseArr(i,j);
    end
  end
end


function [DV] = Dropoff(z)
  %Use this to attenuate the gaussian in z
    DV = -exp(2*const.alpha*(const.z0-z));
end

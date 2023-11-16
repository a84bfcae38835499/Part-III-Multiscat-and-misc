%this function is used to create the potential files.
clear all; close all;
z=linspace(-2,6,100);
%z points, note that the points we take should be evenly spaced
a=2.84;%lattice constant
a1=[a, 0]; a2=[0, a];%base vectors in real space
[b1,b2]=Reciprocal(a1,a2);%reciprocal base lattice vectors
gridp=32;%grid points parallel to the surface in one dimension
i1=1:gridp; i2=1:gridp;%point grids in a1 and a2 directions
V=zeros(length(i1),length(i2),length(z));%potential matrix
D=7.63; alpha=1.1; z0=1.0;%the constants used for generating potential
V0=zeros(size(z));
%now assign values to V0
V0=D*exp(2*alpha*(z0-z))-2*D*exp(alpha*(z0-z));
V0=repmat(reshape(V0,1,1,[]),size(V,1),size(V,2),1);
%repeat matrix V0 to prepare to add it to V
V=V+V0;
for ia1=1:length(i1)
for ia2=1:length(i2)
X(ia1,ia2)=(a1(1)*ia1+a2(1)*ia2)./gridp;%the x grid points
Y(ia1,ia2)=(a1(2)*ia1+a2(2)*ia2)./gridp;%the y grid points
end
end
beta=0.10;%the constant used to create V1
V1=-2*beta*D*exp(2*alpha*(z0-z));
Q=zeros(size(X));
for ia1=1:length(i1)
for ia2=1:length(i2)
Q(ia1,ia2)=cos(2*pi*X(ia1,ia2)/a)+cos(2*pi*Y(ia1,ia2)/a);
for indz=1:length(z)
V(ia1,ia2,indz)=V(ia1,ia2,indz)+V1(indz)*Q(ia1,ia2);
%add V1*Q to V
end
end
end
%===

writematrix(V,"V_boyao.csv")

%===
%now we have got the matrix V, it's time to create the files used by
%Multiscat.
potStructArray.V=V;
%this creates a struct called “potStructArray with a field V,
%and this field is defined by the matrix V in your MATLAB’s Workspace.
Multiscat.PreparePotentialFiles(potStructArray);
%this function creates the pot10001.in file, when executing this line,
%make sure that Multiscat.m is in your current folder.
Multiscat.prepareFourierLabels(V);
%this function creates the FourierLabels.in file, which is also essential
%for Multiscat to get the information about potential.
potStructArray.a1=a1; potStructArray.a2=a2;
%make sure you have the base vectors a1 and a2 in your Workspace.
%The orientation of a1 or a2 does not matter, you just need to get
%the length and angle between them right. As an example, for
%LiF (001) surface, a1=[2.84, 0] and a2=[0, 2.84].
%Note that the length unit is angstrom.
potStructArray.zmin=z(1);
potStructArray.zmax=z(end);
potStructArray.zPoints=length(z);
%here we have an array called z that stores the values of z being
%considered. So zmin and zmax are just the first and the last element
%of z, respectively. The length of z is the number of z points used.
confStruct=Multiscat.createConfigStruct(potStructArray);
Multiscat.prepareConfigFile(confStruct);
%the above two lines create Multiscat.conf
%this function calculates the reciprocal lattice vectors b1 and b2
%from the real space basis vectors a1 and a2
function [b1,b2]=Reciprocal(a1,a2)
a3=[0 0 1];
crsprod=a3*(cross([a1 0],[a2 0]))';
%the volume spanned by a1, a2, and a3
b1=2*pi*cross([a2 0],a3)/crsprod;
b1=b1(1:2);
b2=2*pi*cross(a3,[a1 0])/crsprod;
b2=b2(1:2);
end
%% Plot things! wowie!

data = readmatrix("diffrac10001.out",NumHeaderLines=7,FileType="text")
diffractionChannels = data(:,2:3)
diffractionValues = data(:,4)

%Now we need to make the reciprocal lattice and plot them ig or some shit
%idk. since we know the diffraction channels, I suppose it's a simple case
%of:
% - multiply the diffraction channel by the lattice vector
% - make two matricies storing x and y of those
% - plot the data against that! ez!

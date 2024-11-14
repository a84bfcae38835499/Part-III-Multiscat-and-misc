function Vinterp = interpolate_dft(Nxy,Nz,X,Y,Z,DFTsuper,XDFTsuper,YDFTsuper,ZDFT,VDFTsuper)
  VDFTvect = zeros(DFTsuper*DFTsuper*12*12*19,1);
  %Hard coded values; 12 is the xy resolution of the DFT, 19 is the z resolution
  XDFTvect = VDFTvect;
  YDFTvect = VDFTvect;
  ZDFTvect = VDFTvect;
  index = 0;
  for k = 1:19
    z = ZDFT(k);
     for j = 1:12*DFTsuper
      for i = 1:12*DFTsuper
        if(index + 1 ~= 144*DFTsuper*DFTsuper*(k-1)+12*DFTsuper*(j-1)+i)
          error("Indicies mismatch! Make sure all hard coded values are correct")
        end
        index = 144*DFTsuper*DFTsuper*(k-1)+12*DFTsuper*(j-1)+i;
        disp("index = " + num2str(index))
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
        disp("index2 = " + num2str(index2))
        Xvect(index2) = X(i,j);
        Yvect(index2) = Y(i,j);
        Zvect(index2) = z;
      end
    end
  end
  %Vvect = interp3(XDFTvect,YDFTvect,ZDFTvect,VDFTsuper,Xvect,Yvect,Zvect,'linear');
  Vvect = InterpolatedFn(Xvect,Yvect,Zvect);%<- Beware! this step takes absolutely forever
  if(anynan(Vvect))
    error("Nan found! Consider changing zmax and zmin so it doesn't have to interpolate so far")
  end  
  for k = 1:Nz
    for j = 1:Nxy
      for i = 1:Nxy
        Vinterp(i,j,k) = Vvect(Nxy*Nxy*(k-1)+Nxy*(j-1)+i);
      end
    end
  end
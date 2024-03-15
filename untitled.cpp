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

V = 24.9674*(exp(2*1.1029*(3.1411-z))-2*0.4641*exp(1.1029*(3.1411-z))-2*0.1993*exp(2*0.6477*(3.8323-z))) * ((cos(2*pi*((X/(3.16))-((Y/3.16*sqrt(3))-1)))+cos(4*pi*(Y/3.16*sqrt(3)))+cos(2*pi*((X/(3.16))+(Y/3.16*sqrt(3)))))+3/2)/(4.5)
+ 20.1000*(exp(2*1.1500*(3.2200-z))-2*0.9996*exp(1.1500*(3.2200-z))-2*0.0026*exp(2*1.2439*(4.1864-z))) * ((cos(2*pi*(((X/3.16)-1/2)-((Y/3.16*sqrt(3))-1/2)))+cos(4*pi*(Y/3.16*sqrt(3)))+cos(2*pi*((X/(3.16))+(Y/3.16*sqrt(3)))))+3/2)/(4.5)
+ d*(exp(2*gamma*(z2-z))-2*c*exp(gamma*(z2-z))-2*e*exp(2*lambda*(z3-z))) * (1./( 1+exp((r-mu)*ikbT) )/(1./( 1+exp((r-0.45)*10) )))
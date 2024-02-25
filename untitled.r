  ((backgroundDepth1 * exp(2*alpha1*(z01-z))-2*backgroundDepth1*exp(alpha1*(z01-z))) ...
+ (2*0.1*D0*exp(2*beta1*(z11-z)))) ...
* ((cos(2*pi*((X / (3.16))-(Y / 3.16*sqrt(3))))+cos(4*pi*(Y / 3.16*sqrt(3)))+cos(2*pi*((X / (3.16))+(Y / 3.16*sqrt(3))))) + 3/2)/(4.5) ...
+ ((backgroundDepth2 * exp(2*alpha2*(z02-z))-2*backgroundDepth2*exp(alpha2*(z02-z)))) ...
* ((cos(2*pi*(((X/3.16)-1/2)-((Y/3.16*sqrt(3))-1/2)))+cos(4*pi*(Y / 3.16*sqrt(3))) ...
+ cos(2*pi*((X / (3.16))+(Y / 3.16*sqrt(3))))) + 3/2)/(4.5) ...
+ ((backgroundDepth3 * exp(2*alpha3*(z03-z))-2*backgroundDepth3*exp(alpha3*(z03-z))) ...
+ (2*0.1*D3*exp(2*beta3*(z33-z)))) ...
* ((cos(2*pi*((X / (3.16))-((Y/3.16*sqrt(3))-1)))+cos(4*pi*(Y / 3.16*sqrt(3)))+cos(2*pi*((X / (3.16))+(Y / 3.16*sqrt(3))))) + 3/2)/(4.5);

((backgroundDepth1*exp(2*alpha1*(z01-z))-2*backgroundDepth1*exp(alpha1*(z01-z)))+(2*0.1*D0*exp(2*beta1*(z11-z))))*((cos(2*pi*((X/(3.16))-(Y/(3.16*sqrt(3)))))+cos(4*pi*(Y/(3.16*sqrt(3))))+cos(2*pi*((X/(3.16))+(Y/(3.16*sqrt(3))))))+3/2)/(4.5)+((backgroundDepth2*exp(2*alpha2*(z02-z))-2*backgroundDepth2*exp(alpha2*(z02-z))))*((cos(2*pi*(((X/3.16)-1/2)-((Y/(3.16*sqrt(3)))-1/2)))+cos(4*pi*(Y/(3.16*sqrt(3))-1/2))+cos(2*pi*((X/(3.16)-1/2)+(Y/(3.16*sqrt(3))-1/2))))+3/2)/(4.5)+((backgroundDepth3*exp(2*alpha3*(z03-z))-2*backgroundDepth3*exp(alpha3*(z03-z)))+(2*0.1*D3*exp(2*beta3*(z33-z))))*((cos(2*pi*((X/(3.16))-((Y/(3.16*sqrt(3)))-1)))+cos(4*pi*(Y/(3.16*sqrt(3))-1))+cos(2*pi*((X/(3.16))+(Y/(3.16*sqrt(3))-1))))+3/2)/(4.5);


((backgroundDepth1*exp(2*alpha1*(z01-z))-2*backgroundDepth1*exp(alpha1*(z01-z)))+(2*0.1*D0*exp(2*beta1*(z11-z))))
*((cos(2*pi*((X/(3.16))-(Y/(3.16*sqrt(3)))))+cos(4*pi*(Y/(3.16*sqrt(3))))+cos(2*pi*((X/(3.16))+(Y/(3.16*sqrt(3))))))+3/2)/(4.5)
+((backgroundDepth2*exp(2*alpha2*(z02-z))-2*backgroundDepth2*exp(alpha2*(z02-z))))
*((cos(2*pi*(((X/3.16)-1/2)-((Y/(3.16*sqrt(3)))-1/2)))+cos(4*pi*(Y/(3.16*sqrt(3))-1/2))+cos(2*pi*((X/(3.16)-1/2)+(Y/(3.16*sqrt(3))-1/2))))+3/2)/(4.5)
+((backgroundDepth3*exp(2*alpha3*(z03-z))-2*backgroundDepth3*exp(alpha3*(z03-z)))
+(2*0.1*D3*exp(2*beta3*(z33-z))))*((cos(2*pi*((X/(3.16))-((Y/(3.16*sqrt(3)))-1)))
+cos(4*pi*(Y/(3.16*sqrt(3))-1))+cos(2*pi*((X/(3.16))+(Y/(3.16*sqrt(3))-1))))+3/2)/(4.5);

Qhexfunc(X,Y)
        X_n = X / (const.c);
        Y_n = Y / (const.c*sqrt(3));
        Q = ((cos(2*pi*((X/(const.c))-(Y/(const.c*sqrt(3)))))+cos(4*pi*(Y/(const.c*sqrt(3))))+cos(2*pi*((X/(const.c))+(Y/(const.c*sqrt(3)))))) + 3/2)/(4.5);

          function [V0] = V0func(z,z0,backgroundDepth,alpha)
        V0 = backgroundDepth * exp(2*alpha*(z0-z))...
            -2*backgroundDepth*exp(alpha*(z0-z));
    end
    function [V1] = V1func(z,z0,D,alpha)
        V1 = 2*const.beta*D*exp(2*alpha*(z0-z));


backgroundDepth1*(exp(2*alpha1*(z01-z))-gamma1*exp(alpha1*(z01-z))-2*b*exp(2*beta1*(z11-z))) - A*exp(lambda*-z)
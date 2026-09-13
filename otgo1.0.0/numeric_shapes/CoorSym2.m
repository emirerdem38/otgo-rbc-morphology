function rho = CoorSym2(shp, numDigits, varargin)
% COORSYM2 - Returns the symbolic expressions of rho which is the radial distance
% in polar coordinates.

syms u v;

if ~isempty(varargin)
    radius = varargin{1};
end

switch shp
   case 'dumbbell'
    rho = 1 + real(YnmSym(2,0));

   case 'tiltDumbbell'	
    rho = 1 + real(YnmSym(2,1)) + .1*real(YnmSym(3,2));

   case 'fourBump'
    rho = 1 + .2*exp(-3*real(YnmSym(3,2)));

   case 'eightBump'
    rho = exp(.5*sin(u).^4.*cos(u).*cos(4*v));

   case 'neck'
    rho = 0.1 + real(YnmSym(1,0)).^2;

   case 'sphere'
    rho = radius;

   case 'RBC'

    t_min = 0.81;
    t_max = 2.52;
    r = 3.91;
    d = 2.76;
    
    C0 = t_min ./ (2*r);
    C1 = ((4*r.^2)./(2*d.^2)) .* (-4*C0 + ((t_max.*abs(5*d.^2-16.*r.^2))./((4*r.^2-d.^2).^(3/2))));
    C2 = ((16*r.^2)./(2*d.^4)) .* (2*C0 + ((t_max.*(16*r.^2-3*d.^2))./(sign(5*d.^2-16*r.^2).*(16*r.^2-d.^2).^(3/2))));
    
    rho = ((1-sin(u).^2).*(C0+C1*sin(u).^2+C2*sin(u).^4)).^(1/2.*sign(cos(u)));
      
end

rho = vpa(rho, numDigits);

end


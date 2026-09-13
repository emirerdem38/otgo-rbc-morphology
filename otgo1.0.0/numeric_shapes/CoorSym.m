function [X, Y, Z] = CoorSym(shp)
% COORSYM - Returns the symbolic expressions of x,y,z components of a cell with
%  the shape SHP.

syms u v
ax = 1; ay = 1; az = 1;

switch shp
   case 'ellipseZ'
    rho = 1;
    ax = 1/2; ay = 1/2;

   case 'ellipseX'
    rho = 1;
    ay = .34;
    az = .34;

   case 'dumbbell'
    rho = 1 + real(YnmSym(2,0));

   case 'tiltDumbbell'	
    rho = 1 + real(YnmSym(2,1)) + .1*real(YnmSym(3,2));

   case 'fourBump'
    rho = 1 + .2*exp(-3*real(YnmSym(3,2)));

   case 'eightBump'
    rho = exp(.5*sin(u).^4.*cos(u).*cos(4*v));

   case 'oblate85'
    rho = 1;
    ax = .47;

   case 'oblate75'
    rho = 1;
    ax = .366;

   case 'oblate65'
    rho = 1;
    ax = .29;
    
   otherwise
    rho = 1;
end

X = ax*rho*sin(u)*cos(v);
Y = ay*rho*sin(u)*sin(v);
Z = az*rho*cos(u);

end
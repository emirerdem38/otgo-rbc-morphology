function YnSym = YnmSym(n,m)
% YNSYM - Returns the symbolic expression for spherical harmonics of degree
% n and order m.
%
    PnSym = [];
    syms v
    for i = 0:n
        PnSym = [PnSym, LegendreSym(n, i)];
    end
    
    nor = getNorConst(n);
    PnSym = nor.*PnSym(:,[n+1:-1:2 1:n+1]);
    YnSym = PnSym.*exp(1i*v*(-n:n));
    YnSym = YnSym(:,n+1+m);
end

function nc = getNorConst(n)
% GETNORCONST - returns the normalizing constant for the Legendre
% polynomials of degree n
%
  persistent NOR

  if(n>=length(NOR) || isempty(NOR) || isempty(NOR{n+1}))
    NOR{n+1} = [(-1).^(-n:0) ones(1,n)].*sqrt((2*n+1)/4/pi *...
                                              myFact(n-abs(-n:n))./myFact(n+abs(-n:n)));
  end
  nc = NOR{n+1};
end

function val = myFact(N)
% MYFACT - tabulated factorial
%
  persistent factTable

  q = max(N);
  if(q>=length(factTable))
    tail = length(factTable);
    factTable = [factTable factorial(tail:q)];
  end
  val = factTable(N+1);
end













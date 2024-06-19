F<w>:=QuadraticField(13);
OF:=Integers(F);

I2:=2*OF;
I13:=Factorisation(13*OF)[1,1];
I3a:=Factorisation(3*OF)[1,1];
I3b:=Factorisation(3*OF)[2,1];
R<e>:=PolynomialRing(Rationals());


// primes of F
PRIMEBOUND:=200;

load "PrimeIdeals13.out";

PrimeIdeals:=[];
for pi in PrimeGenerators do
  Append(~PrimeIdeals,(OK!pi)*OK);
end for;


procedure PrintNewforms(filename,label,forms);

PrintFile(filename, label cat ":=[":Overwrite:=true);

for i:=1 to #forms do
  f:=forms[i];

  EI:=[];
  for Q in PrimeIdeals do
    Append(~EI,HeckeEigenvalue(f,Q));
  end for;

  KF:=Parent(EI[1]);
  FI:=<R!DefiningPolynomial(KF),EI>;
  
  PrintFile(filename,FI);
  if i lt #forms then
    PrintFile(filename,",");
  end if;
end for;

PrintFile(filename, "];");

end procedure;


/*
N1:=I2^3*I13^2;
print "computing newforms of level N1. Space of dimension ", Dimension(NewSubspace(HilbertCuspForms(K,N1)));
time NewformsN1:=Eigenforms(NewSubspace(HilbertCuspForms(K,N1)));
print "..done.";
PrintNewforms("Newforms13-E-N1.out","CoefficientsN1",NewformsN1);
*/


N2:=I2^4*I13^2;
print "computing newforms of level N2. Space of dimension ", Dimension(NewSubspace(HilbertCuspForms(F,N2)));
time NewformsN2:=Eigenforms(NewSubspace(HilbertCuspForms(F,N2)));
print "..done.";
PrintNewforms("Newforms13-E-N2.out","CoefficientsN2",NewformsN2);







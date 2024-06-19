
K<w>:=QuadraticField(13);
OK:=Integers(K);
P<x>:=PolynomialRing(Integers());
R<e>:=PolynomialRing(Rationals());
 
// save common list of primes of K
PRIMEBOUND:=200;
PrimeIdeals:=[Q: Q in PrimesUpTo(PRIMEBOUND,K) | Norm(Q) mod 2 ne 0 and Norm(Q) mod 13 ne 0];

filename:="PrimeIdeals13.out";

PrintFile(filename,"PrimeGenerators:=[": Overwrite:=true);

for i in [1..#PrimeIdeals] do
  Q:=PrimeIdeals[i];
  _,g:=IsPrincipal(Q);
  PrintFile(filename,K!g);
  if i lt #PrimeIdeals then
    PrintFile(filename,",");
  end if;
end for;

PrintFile(filename,"];");


/*

In this file, we compute a list of prime ideals of Q(sqrt(13)) of norm up to 200.
We store them in the file 13-curveE_PrimeIdeals.out.

*/


// The field F = Q(sqrt(13)) (i.e., the unique quadratic subfield in Q(zeta_13))
F<w>:=QuadraticField(13);
OF:=Integers(F);
P<x>:=PolynomialRing(Integers());
R<e>:=PolynomialRing(Rationals());
 
// save common list of primes of F = Q(sqrt(13))
PRIMEBOUND:=200;
PrimeIdeals:=[Q: Q in PrimesUpTo(PRIMEBOUND,F) | Norm(Q) mod 2 ne 0 and Norm(Q) mod 13 ne 0];

filename:="13-curveE_PrimeIdeals.out";

PrintFile(filename,"PrimeGenerators:=[": Overwrite:=true);

for i in [1..#PrimeIdeals] do
    Q:=PrimeIdeals[i];
    _,g:=IsPrincipal(Q);
    PrintFile(filename,F!g);
    if i lt #PrimeIdeals then
        PrintFile(filename,",");
    end if;
end for;

PrintFile(filename,"];");


/*

This file contains the computations to prove Proposition 9 for the level 2^4w^2.

*/

load "13-curveE.m";


N42:=I2^4*I13^2;

print "Computing newforms of level 2^4*w^2. Space of dimension", Dimension(NewSubspace(HilbertCuspForms(F,N42)));
time forms42:=Eigenforms(NewSubspace(HilbertCuspForms(F,N42)));
print "...done!";
print "There are",#forms42,"newforms to eliminate.";


BoundE(forms42, [3, 17, 23, 29, 43]);


print "We now check the remaining form at the level 2^4*w^2 (i.e. when s=4 in Lemma 7) corresponds to E(1,1) by comparing Hecke coefficients at primes of norm up to 200:";
W0:=FreyE(1,1);
print {HeckeEigenvalue(forms42[67],Q) - TraceOfFrobenius(W0,Q) : Q in PrimesUpTo(200,F) | Valuation(N42,Q) eq 0} eq {0};
print "Done!";
print "This finishes the proof of Proposition 9 in the large level case.";


  

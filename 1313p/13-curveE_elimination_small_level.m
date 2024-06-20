/*

This file contains the computations to prove Proposition 9 at level 2^3*w^2.

*/

load "13-curveE.m";

N32:=I2^3*I13^2;

print "Computing newforms of level 2^3*w^2. Space of dimension", Dimension(NewSubspace(HilbertCuspForms(F,N32)));
time forms32:=Eigenforms(NewSubspace(HilbertCuspForms(F,N32)));
print "...done!";
print "There are",#forms32,"newforms to eliminate.";



BoundE(forms32, [3,17,23,29]);

W0:=FreyE(1,-1);
W1:=FreyE(1,0);

print "We check that the surviving forms (for all p different from 7) at the level 2^3*w^2 (i.e. when s = 3 in Lemma 7) correspond to solutions (1,-1,0) and (1,0,1).";
assert {HeckeEigenvalue(forms32[12],Q) - TraceOfFrobenius(W0,Q) : Q in PrimesUpTo(200,F) | Valuation(N32,Q) eq 0} eq {0};
assert {HeckeEigenvalue(forms32[13],Q) - TraceOfFrobenius(W1,Q) : Q in PrimesUpTo(200,F) | Valuation(N32,Q) eq 0} eq {0};
print "Done!";

print "";
print "For p = 7 the surviving forms correspond to the trivial solutions (1,-1,0) and (1,0,1) plus one extra form with quadratic field of coefficients Q(Sqrt(2)).";
print "This completes the proof of Proposition 9 for the level 2^3*w^2.";
print "";
print "";

print "We now check Remark 7.4.";
print "Let g be the form not eliminated for p = 7";
f27:=forms32[27];
Qf<s>:=BaseField(f27);
assert IsIsomorphic(Qf,QuadraticField(2));
Of:=Integers(Qf);
p7,p7prime:=Explode([p[1] : p in Factorisation(7*Of)]);
assert (s + 3) in p7;
assert (s + 4) in p7prime;

print "";
print "Do the traces at Q of g and E(1,-1) match mod P7' for Q of norm up to 5000?";
k,tok:=ResidueClassField(p7prime);
{tok(HeckeEigenvalue(f27,Q) - TraceOfFrobenius(W0,Q)) : Q in PrimesUpTo(200,F) | Valuation(N32,Q) eq 0} eq {0};

print "";
print "Do the traces at Q of g and E(1,-1) match mod P7 for Q of norm up to 5000?";
k,tok:=ResidueClassField(p7);
{tok(HeckeEigenvalue(f27,Q) - TraceOfFrobenius(W0,Q)) : Q in PrimesUpTo(200,F) | Valuation(N32,Q) eq 0} eq {0};

print "";
print "This non rational form seems to have a mod 7 representation isomorphic to that of the Frey curve attached to the trivial solution (1,-1,0)";



  

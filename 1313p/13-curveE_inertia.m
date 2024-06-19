/*

This file contains the computations to prove (part (B) of) Theorem 7 which is stated on page 8665 and proved on page 8668.

*/


load "13-curveE.m";



// Definition of the field M

Z:=FreyE(1,-1);
f3a:=Factorization(DivisionPolynomial(Z,3))[1][1];
F2<x>:=NumberField(f3a);
R<y>:=PolynomialRing(F2);
M:=ext<F2|Evaluate(DefiningPolynomial(Z),[x,y,1])>;

print "Choice of 3-torsion point whose coordinates (x,y) satisfy:";
print f3a;
print "Adjoining y-coordinate:", Evaluate(DefiningPolynomial(Z),[x,y,1]);    
print "Degree of M over F =", AbsoluteDegree(M)/Degree(F);


OM:=Integers(M);

J:=Factorization(2*OM);

assert #J eq 1;
print "";
print "There is a unique prime above 2 in M.";

J2:=J[1][1];

print "";
print "It is totally ramified (with ramification index 8). We call it Q'.";



// Inertia argument over M.

S:={};
print "";
print "We display below a, b (within {1,..8}, not both even and such that 4 does not divide a + b), the valuation of the conductor of E(a,b) base changed to M at Q' and the Kodaira Symbol.";
for a,b in [1..2^3] do
    if (a le b) and ([a mod 2, b mod 2] ne [0,0]) and ((a+b) mod 4 ne 0) then
    E1:=FreyE(a,b);
    LI:=LocalInformation(BaseChange(E1,M),J2);
    L1:=LI[3];
    print "a =",a,"b =",b,"val =",L1, "Kodaira type", LI[5];
    S:=S join {L1};
    end if;   
end for;

print "";
print "The valuation of the conductor of E(a,b) at Q' belongs to:";
print S;





// test if fixed elliptic curve over M has conductor exponent 2 at a prime over 2
print "We check that the conductor exponent of Z = E(1,-1) base changed to M at the unique prime ideal Q' of M above 2 is 2.";
assert LocalInformation(BaseChange(Z,M),J2)[3] eq 2;
print "Done!";


print "This proves Theorem 7.";



print "In addition, we check that the full 3-division field of Z = E(1,-1) has degree 48 over Q(sqrt{13}).";

R<t>:=PolynomialRing(F);

j:=jInvariant(Z);
// Has image in normalizer of non-split Cartan?
assert Roots(j - t^3) eq [];

// Has image in Borel?
assert Roots(j*t^3 - (t+27)*(t+243)^3) eq [];

// Has image in normalizer of split Cartan?
assert Roots(j*t^3 - ((t-9)*(t+3))^3) eq [];

print "Done!";
/*
We prove that the mod p reperesentation attached to F is irreducible for p > 5 and p different from 13 when either 13 |  a + b 
We prove that the mod 5 representation attached to F is irreducible. 
See Lemma 10 p. 8672 in the published version of the paper.

*/


load "13-curveF.m";



print "We compute the integer B that appears in Theorem 1 in Freitas-Siksek's 2015 paper (JTNB).";
print "";

S:= {[0,0,12], [0,12,0], [12,0,0], [0,12,12], [12,0,12], [12,12,0]}; // list of possible signatures

// The unit group in K (= unique cubic subfield in Q(zeta13)))
U,phi:=UnitGroup(K);

G:=Automorphisms(K); // Automorphism group of K

// This function computes the twisted norm associated to a given signature as in Freitas-Siksek's paper
function TwistedNorm(a,s);
    N:=1;
    i:=1;
    for g in G do
        N:=N*g(a)^s[i];
        i:= i + 1;
    end for;
    return N;
end function;

// This function computes the number A_s that appears in Theorem 1 of Freitas-Siksek
function As(s);
    return Norm(Gcd(Gcd((TwistedNorm(K!phi(U.1),s) - 1)*OK,(TwistedNorm(K!phi(U.2),s) - 1)*OK),(TwistedNorm(K!phi(U.3),s) - 1)*OK));
end function;


B:=Lcm([As(s) : s in S]);

print "The integer B for the number field K from Theorem 1 in Freitas-Siksek (JTNB, 2015) is:";
print B;
print "Its prime factors are:";
print PrimeFactors(B);


I2:=Factorisation(2*OK)[1,1];
I13:=Factorisation(13*OK)[1,1];


print "";
print "The Ray Class Group for modulus Q2^2*oo1*oo2*oo3 is:";
R,_:=RayClassGroup(I2^2,[1,2,3]);
print R;


print "";
print "The Ray Class Group for modulus Q13*oo1*oo2*oo3 is:";
R,_:=RayClassGroup(I13,[1,2,3]);
print R;


print "";
print "The Ray Class Group for modulus Q2^2*Q13*oo1*oo2*oo3 is:";
R,_:=RayClassGroup(I2^2*I13,[1,2,3]);
print R;

// Definition of the Frey curve F as

L<x,y>:=PolynomialRing(K,2);
KF<a,b>:=FieldOfFractions(L);


F:=FreyF(a,b);

print("Numerator of j(F) - 1728:");
Numerator(jInvariant(F)-1728);

print("Factorisation of the numerator of j(F) - 1728:");
Factorization(Numerator(jInvariant(F)-1728));



print("Denominator of j(F) - 1728:");
Denominator(jInvariant(F)-1728);

print("Factorisation of the denominator of j(F) - 1728:");
Factorization(Denominator(jInvariant(F)-1728));



print("Definition of eta:");
eta:=4864*z^2 - 7168*z + 1600;
print(eta);


print("Is eta/13 a square?");
IsSquare(eta/13);

// Definition of D : 13y^2 = (x^2 + 22*x + 125)*x *over K*  
SK<x>:=PolynomialRing(K);
DK:=QuadraticTwist(EllipticCurve((x^2 + 22*x + 125)*x),13);

//DK has rank 1 over K
print("Rank of D/K");
RankBounds(DK);

//DK has Z/2 torsion subgroup
print("Torsion subgroup of D/K");
TorsionSubgroup(DK);

//DK in fact has a model DQ over Q.
//DQ in fact has rank 1 over Q.
//DQ has Z/2 torsion subgroup so torsion is already defined over Q.

// Definition of D : 13y^2 = (x^2 + 22*x + 125)*x *over Q*  
DQ:=QuadraticTwist(EllipticCurve((X^2 + 22*X + 125)*X),13);
print("Rank of D/Q");
RankBounds(DQ);





// Check when j(F) - 1728 is rational.

G:=Automorphisms(K);

nj:=Evaluate(Numerator(jInvariant(F)-1728),[x,1]);
dj:=Evaluate(Denominator(jInvariant(F)-1728),[x,1]);

function apply_galois(g,f)
    c1:=Coefficients(f);
    c2:=[g(x) : x in c1];
    return elt<SK| c2>;
end function;

NJ:=nj*apply_galois(G[2],dj)*apply_galois(G[3],dj);

print("Here is R:");
print(NJ);

function extract_equations(f)
    e1:=0;
    e2:=0;
    e3:=0;
    for i:=0 to Degree(f) do
        c:=Coefficient(f,i);
        e1:=e1 + x^i*c[1];
        e2:=e2 + x^i*c[2];
        e3:=e3 + x^i*c[3];
    end for;
    return SK!e1,SK!e2,SK!e3;
end function;


e1,e2,e3:=extract_equations(NJ);



(e1 + e2*z + e3*z^2) eq NJ;

print "The polynomial A1 is";
e2;
print "Its roots are:";
Roots(e2);

print "The polynomial A2 is";
e3;
print "Its roots are:";
Roots(e3);

print "The only common root is x = -1, which corresponds to (a,b) = +/-(1,-1).";

print "This proves Lemma 10.";


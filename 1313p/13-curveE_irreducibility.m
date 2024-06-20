/*

We prove that the mod p representation attached to E/Q(sqrt(13)) is irreducible for p > 5
and that the mod 5 representation is also irreducible when 3 | a + b.
See Proposition 8 p. 8667 in the published version of the paper.

*/

load "13-curveE.m";

P<x>:=PolynomialRing(Integers());

print "We apply Theorem 1 in Freitas-Siksek's 2015 paper (JTNB) using q = 3 and 17.";
print "";

S:= {[0,12], [12,0]}; // list of possible signatures

// The unit group in F
U,phi:=UnitGroup(F);

G:=Automorphisms(F); // Automorphism group of F

// This function computes the twisted norm associated to a given signatures as in Freitas-Siksek's paper
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
    return Norm(Gcd((TwistedNorm(F!phi(U.1),s) - 1)*OF,(TwistedNorm(F!phi(U.2),s) - 1)*OF));
end function;


B:=Lcm([As(s) : s in S]);

print "The integer B for the number field Q(sqrt(13)) from Theorem 1 in Freitas-Siksek (JTNB, 2015) is:";
print B;
print "Its prime factors are:";
print PrimeFactors(B);

print "";
print "We now apply Theorem 1 in Freitas-Siksek (JTNB, 2015) using q = 3, 17 (as in the proof of Theorem 3 in loc. cit.).";
print "";

factors:=[];
for q in [3,17] do
    q1:=Factorisation(q*OF)[1,1];
    q2:=Factorisation(q*OF)[2,1];

    tracesq:=[];

    for a,b in [0..(q-1)] do 
        if (a le b) and ([a,b] ne [0,0]) then
            C:=FreyE(a,b); 
            locQ1:=LocalInformation(C,q1);
            locQ2:=LocalInformation(C,q2);
            assert locQ1[3] eq 0 and locQ2[3] eq 0; // The curve E has good reduction at both primes above q
            Append(~tracesq,[TraceOfFrobenius(C,q1),TraceOfFrobenius(C,q2)]);
        end if;
    end for;

    Append(~factors,PrimeDivisors(Lcm([Gcd([Resultant(x^2 - tr[i]*x + q, x^12 - 1) : i in [1,2]]) : tr in tracesq])));
end for;

//assert #factors eq 2;
print "The set of primes p such that p = 11 or p > 13 (and p does not divide B) for which reducibility may occur is:";
S:=Set(factors[1]) meet Set(factors[2]);
S:={p : p in S | p notin {2,3,5,7,13} and B mod p ne 0};
S;
print "Using Theorem 1 in Freitas-Siksek (JTNB, 2015), this proves that irreducibility is missing only for p = 5,7,13.";


print "";
print "We already know that E has good reduction at both primes Q1 and Q2 above 3 in Q(sqrt(13)). We compute the traces:";
Q1:=Factorisation(3*OF)[1,1];
Q2:=Factorisation(3*OF)[2,1];

traces3:=[];
for a,b in [0..2] do 
  if (a le b) and ([a,b] ne [0,0]) then
	C:=FreyE(a,b); 
        locQ1:=LocalInformation(C,Q1);
        locQ2:=LocalInformation(C,Q2);
        print "If a, b (mod 3) =",a,",",b,"then (aQ1(E(a,b), aQ2(E(a,b))) = ",[TraceOfFrobenius(C,Q1),TraceOfFrobenius(C,Q2)];
	assert locQ1[3] eq 0 and locQ2[3] eq 0;
        Append(~traces3,[TraceOfFrobenius(C,Q1),TraceOfFrobenius(C,Q2)]);
   end if;
end for;

traces3:=SetToSequence(Set(traces3));

print "For each possible pair of traces (aQ1(E(a,b), aQ2(E(a,b))) we check whether the polynomials";
print "(X^2 - aQ1(E(a,b))X + 3 (mod p), X^2 - aQ1(E(a,b))X + 3 (mod p)) for p = 5, 7, 13";
print "are irreducible.";

for p in [5,7,13] do
    Py<y>:=PolynomialRing(GF(p));
    for tr in traces3 do
        print "Dealing with p =",p,"and traces =",tr;
        [IsIrreducible(y^2 - tr[i]*y + 3) : i in [1,2]];
    end for;
end for;

print "This shows that the mod p representation attached to E is also irreducible for p = 5,7,13 except in the case p = 5 and trace of Frobenius at 3 equal to [-1,1].";
print "Since if 3 | a + b, we have traces [-3,-1]. The result follows.";

print "This proves Proposition 8.";



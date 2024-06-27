/*

This file contains the computations to prove Proposition 9.

*/

load "../13-curveE.m";

R<e>:=PolynomialRing(Rationals());


// loads the list "PrimeGenerators" which contains generators for primes in K with norm absolute value up to 200
// the primes above 2 and 13 are not included, as they appear in the level of the modular forms and we are interested in taking traces of Frobenius
load "13-curveE_PrimeIdeals.out";

// creates primes of F
PrimeIdeals:=[];
for pi in PrimeGenerators do
  Append(~PrimeIdeals,(OF!pi)*OF);
end for;

// the space of forms takes a while to compute, so we have precomputed the coefficients of each form for the primes in PimeIdeals (in that order)
// load coefficients of forms of level I2^4*I13^2 -- this is case s = 4 in Lemma 7
load "13-curveE_Newforms_large_level.out";



// Returns the set of prime divisors of x, with special treatment when x is zero or a unit.
S:=Parent({1,2,3});
T:=Parent({{1,2,3}});

function primeset(x)
    x:=Integers()!x;
    if x eq 0 then
        return S!{0};
    elif x in {-1,1} then
        return S!{x};
    else 
        return S!Set(PrimeDivisors(x));
    end if;
end function;

// Find the index in PrimeIdeals that corresponds to the prime ideal Q.
function LookUpIndex(Q);
    for i:=1 to #PrimeIdeals do
        if Q eq PrimeIdeals[i] then
            return i;
        end if;
    end for;
end function;

// returns the field of coefficients of the form cf, together with a primitive element s
function GetNumberField(cf);
    if Degree(cf[1]) eq 1 then
        Kf:=Rationals();
        s:=-Coefficient(cf[1],0);
    else
        Kf<s>:=NumberField(cf[1]);
    end if;
    return Kf,s;
end function;


// returns true if form cf and elliptic curve W have the same traces at all primes in PrimeIdeals, otherwise returns false
function SameTraces(cf,W);
    Kf,s:=GetNumberField(cf);
    return Set([Evaluate(cf[2][LookUpIndex(Q)],s) - TraceOfFrobenius(W,Q) : Q in PrimeIdeals]) eq {0};
end function;

// Let cf be a modular form with a quadratic field of coefficients Kf and p a prime splitting in Kf
// tests if the form cf and elliptic curve W have the same traces mod P (at all primes in PrimeIdeals coprime to p) for the primes dividing p

function SameTracesModP(cf,W,p);
    Kf,s:=GetNumberField(cf);
    assert Degree(Kf) eq 2;	

    OFf:=Integers(Kf);

    // Pick a prime above p.
    factp:=Factorization(p*OFf);	 
    assert #factp eq 2;	
    P1:=factp[1][1];
    P2:=factp[2][1];
    k1,pi1:=ResidueClassField(P1);
    k2,pi2:=ResidueClassField(P2);

    A1:=Set([pi1(Evaluate(cf[2][LookUpIndex(Q)],s) - TraceOfFrobenius(W,Q)) : Q in PrimeIdeals | Norm(Q) mod p ne 0]) eq {k1!0};
    A2:=Set([pi2(Evaluate(cf[2][LookUpIndex(Q)],s) - TraceOfFrobenius(W,Q)) : Q in PrimeIdeals | Norm(Q) mod p ne 0]) eq {k2!0};

    return P1,A1,P2,A2;
end function;

// This function computes the bounds \mathcal{C}_q(f) for C = E
// it assumes that the curve has good reduction at the primes above q
function Bound(q,cf,curve)
    B:=1;

    Kf,s:=GetNumberField(cf);

    factQ:=Factorisation(q*OF);
    for x,y in [0..q-1] do
        if [x,y] ne [0,0] then
            L:=[];
            C:=curve(x,y);
            for i in [1..#factQ] do
                Q:=factQ[i,1];
                afQ:=Evaluate(cf[2][LookUpIndex(Q)],s);
                L:=Append(L,Integers()!Norm(TraceOfFrobenius(C,Q)-afQ));
            end for;
            B:=B*Gcd(L);
        end if;
    end for;
    return primeset(q*B);
end function;


// Eliminate forms at level I2^4*I13^2 and compute the bounds on p.
// note the auxiliary primes used are all not congruent to 1 mod 13, so of good reduction for the Frey curve by Lemma 5

cforms:=CoefficientsN2; 
curve:=FreyE;
survE:=[];
Boundsall:={};
print "We eliminate the forms at larger level 2^4*w^2 (s = 4). There are", #cforms, "forms to eliminate.";
for i in [1..#cforms] do
    print "************** Dealing with form no",i,"*************";
    cf:=cforms[i];
    Bf:={0};
    bool:=0;
    for q in [3,17,23,29,43,61] do
        if (bool eq 0) or (Bf ne {}) then
            print "Dealing with auxiliary prime q =",q;
            Bqf:=Bound(q,cf,curve);
            //print Bqf;
            Bqf:={x : x in Bqf | x notin {2,3,13}};
            if Bqf ne {0} then
                if bool eq 0 then
                    print "This form can be eliminated for large enough p !";
                    bool:=1;
                    Bf:=Bqf;
                end if;
                Bf:= Bf meet Bqf;
                print "Prime exponent(s) remaining to eliminate:",Bf;
            end if;
            //Bf:=Append(Bf,Bqf);
        end if;
    end for;
    
    if Bf eq {0} then
        Append(~survE,i); 
        print "Form no",i,"not eliminated!";
    else
        //Boundsf:=&meet {x : x in Bf | x ne {0}};
        //Boundsf:={x : x in Boundsf | x notin {2,3,13}};
        if #Bf eq 0 then
            print "Form no",i,"eliminated!";
        else
            print "Form no",i,"eliminated for exponents not in",Bf;
        end if;
        //Boundsall:=Boundsall join Boundsf;
    end if;
end for;
//assert Boundsall eq {};
print "Done!";
print "+++++++++++++++++++++++++++++++++++++++";


assert #survE eq 1;
W:=FreyE(1,1);
assert SameTraces(cforms[survE[1]],W);
print "The surviving form at level 2^4*w^2 (s = 4) corresponds to (1,1).";
print "+++++++++++++++++++++++++++++++++++++++";
print "This proves Proposition 9 in the 'large' level 2^4*w^2.";

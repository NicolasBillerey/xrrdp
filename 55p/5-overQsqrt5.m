/*

We define here the field K = Q(sqrt(5)) and all the useful functions for elimination.
We also compute the weight-2 newforms of level 2^6 (which are used for Proposition 4) over K and define the Frey curves E/K and F/K.

*/

K<s>:=QuadraticField(5);
OK:=Integers(K);
o:=(-1 + s)/2;
ob:=(-1 - s)/2;

I2:=Factorisation(2*OK)[1,1];
I3:=Factorisation(3*OK)[1,1];
//I5:=Factorisation(5*OK)[1,1];

// This space of forms is used in the proof of Proposition 4.
N:=I2^6;
print "Computing newforms over Q(sqrt(5)) of weight 2 and level 2^6...";
Nfs:=Eigenforms(NewSubspace(HilbertCuspForms(K,N)));
print "...done!";


// The following is used in the proof of Proposition 4.

PrimeIdeals:=[Q: Q in PrimesUpTo(30,K) | Norm(Q) mod 2 ne 0];

function SameTraces(f,curve);
  return Set([HeckeEigenvalue(f,Q) - TraceOfFrobenius(curve,Q) : Q in PrimeIdeals]) eq {0};
end function;

// Definition of the elliptic curve E/Q(sqrt5) and its twists
function FreyE(a,b,d);
    return EllipticCurve([0,2*(a+b)*d,0,-ob*(a^2 + o*a*b + b^2)*d^2,0]);
end function;

// Definition of the elliptic curve F/Q(sqrt5) and its twists
function FreyF(a,b,d);
    return EllipticCurve([0,2*(a-b)*d,0,(-3*s/10 + 1/2)*(a^2 + o*a*b + b^2)*d^2,0]);
end function;


// This function is used to eliminate (exponents for) f using the prime q and NORMS of differences of traces
function BoundForm(q,f,curve,exponents,bound)
    // q is the auxiliary prime
	// f is the form to eliminate
	// exponents should be a set with the prime exponents to eliminate; if no restrictions are known set exponents:={}
    // bound is a lower bound on the exponents to eliminate, that is only prime exponents p greater than or equal to the bound will be considered.
    // Note that the ouput is {} iff q does NOT bound any exponent.
    // Moreover, if the ouput is different from {}, then it contains q by construction, but this is harmless for the elimination procedure.
    B:=1;
    factQ:=Factorisation(q*OK);
    for x,y in [0..q-1] do
        if ([x,y] ne [0,0]) and (x le y) then
            Bxy:=0; 
            C:=curve(x,y,1);
            for i in [1..#factQ] do
                Q:=factQ[i,1];
                if LocalInformation(C,Q)[3] eq 0 then
                    // Here C has good reduction at Q
                    diffQ:=Norm(TraceOfFrobenius(C,Q) - HeckeEigenvalue(f,Q));
                else
                    diffQ:=Norm((Norm(Q)+1)^2 - HeckeEigenvalue(f,Q)^2);
                end if;
                Bxy:=Gcd(Bxy,Integers()!diffQ);
            end for;
            if Bxy eq 0 then
            	return {}; // Here p is unbounded
            else
              B:=B*Bxy;
            end if;
        end if;
    end for;
	assert B ne 0;
    Sf:={p : p in Set(PrimeDivisors(B)) | p ge bound};
    if exponents ne {} then
      return (Sf meet exponents) join {q};
    else
      return Sf join {q};
    end if;
end function;


// This function is used to eliminate each form in a given space using specified auxiliary primes and NORMS of differences of traces
function Bound(forms,curve,AuxiliaryPrimes,bound);
	print "Performing standard elimination for",#forms,"form(s) with set of auxiliary primes",AuxiliaryPrimes;
	for i in [1..#forms] do
		f:=forms[i];
		print "";
		print "Checking form no",i;
		print "";
		Sf:={};
		bool:=0;
        for q in AuxiliaryPrimes do
            if bool eq 0 or Sf ne {} then
                print "Dealing with q =",q;
                Sq:=BoundForm(q,f,curve,Sf,bound);
                if Sq ne {} then // Here f can be discarded for large enough p
                    if bool eq 0 then
                        print "This form can be eliminated for large enough p !";
                        Sf:=Sq;
                        bool:=1;
                    end if;
                    Sf:=Sf meet Sq;
                    Sf:={p : p in Sf | p ge bound};
                end if;
                if Sf ne {} then
                    print "Prime exponents remaining to eliminate =",Sf;
                end if;
            end if;
        end for;
		if bool eq 0 then
		    print "Form no",i,"not eliminated for large enough p";
        else
            if Sf eq {} then
                print "Form no",i,"is eliminated";
            else
                print "Form no",i;
                print "with coefficient field :", CoefficientField(f) ;
                print "is not eliminated for prime(s) :",Sf;
            end if;
		end if;
		print "*************************************************************";
	end for;
	return "";
end function;



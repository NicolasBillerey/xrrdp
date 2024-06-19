/*

This file contains the computations to prove Proposition 9 for the level 2^4w^2.

*/

load "13-curveE.m";


I2:=Factorization(2*OF)[1,1]; 
I13:=Factorization(13*OF)[1,1];

N42:=2^4*I13^2;

print "Computing newforms of level 2^4*w^2. Space of dimension", Dimension(NewSubspace(HilbertCuspForms(F,N42)));
time forms42:=Eigenforms(NewSubspace(HilbertCuspForms(F,N42)));
print "...done!";
print "There are",#forms42,"newforms to eliminate.";


// This function is used to eliminate (exponents for) f using the prime q and NORMS of differences of traces
function BoundFormE(q,f,exponents)
    // q is the auxiliary prime
	// f is the form to eliminate
	// exponents should be a set with the prime exponents to eliminate; if no restrictions are known set exponents:={}
    // Note that the ouput is {} iff q does NOT bound any exponent.
    // Moreover, if the ouput is different from {}, then it contains q by construction, but this is harmless for the elimination procedure.
    B:=1;
    factQ:=Factorisation(q*OF);
    for x,y in [0..q-1] do
        if ([x,y] ne [0,0]) and (x le y) then
            Bxy:=0; 
            C:=FreyE(x,y);
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
    Sf:={p : p in Set(PrimeDivisors(B)) | p notin {2,3,13}};
    if exponents ne {} then
      return (Sf meet exponents) join {q};
    else
      return Sf join {q};
    end if;
end function;


// This function is used to eliminate each form in a given space using specified auxiliary primes and NORMS of differences of traces
function BoundE(forms,AuxiliaryPrimes);
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
                Sq:=BoundFormE(q,f,Sf);
                if Sq ne {} then // Here f can be discarded for large enough p
                    if bool eq 0 then
                        print "This form can be eliminated for large enough p !";
                        Sf:=Sq;
                        bool:=1;
                    end if;
                    Sf:=Sf meet Sq;
                    Sf:={p : p in Sf | p notin {2,3,13}};
                end if;
                print "Prime exponents remaining to eliminate =",Sf;
            end if;
        end for;
		if bool eq 0 then
		    print "Form no",i," is NOT eliminated for large enough p";
        else
            if Sf eq {} then
                print "Form no",i,"is eliminated";
            else
                print "Form no",i;
                print "with coefficient field :", CoefficientField(f) ;
                print "is NOT eliminated for prime(s) :",Sf;
            end if;
		end if;
		print "*************************************************************";
	end for;
	return "";
end function;

BoundE(forms42, [3, 17, 23, 29, 43, 61]);

/*
To be added when we know which form number corresponds to W0
W0:=FreyE(1,1);
assert {HeckeEigenvalue(forms42[??],Q) - TraceOfFrobenius(W0,Q) : Q in PrimesUpTo(200,F) | Valuation(N42,Q) eq 0} eq {0};
print "The surviving form at the level 2^4*w^2 (i.e. when s=4 in Lemma 7) correspond to the Frey curve evaluated at (1,1).";
*/


  

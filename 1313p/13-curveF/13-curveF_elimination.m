/*

This file contains the computations to prove Theorem 2.

*/


load "13-curveF.m";
load "13-curveE.m";

I2:=Factorization(2*OK)[1,1];
I3:=Factorization(3*OK)[1,1];
I13:=Factorization(13*OK)[1,1];

N111:=I2*I3*I13;

print "Computing newforms of level Q2*Q3*Q13. Space of dimension", Dimension(NewSubspace(HilbertCuspForms(K,N111)));
time forms111:=Eigenforms(NewSubspace(HilbertCuspForms(K,N111)));
print "...done!";
print "There are",#forms111,"newforms to eliminate.";


// This function is used to eliminate (exponents for) f using the prime q and NORMS of differences of traces
function BoundFormF(q,f,exponents)
    // q is the auxiliary prime
	// f is the form to eliminate
	// exponents should be a set with the prime exponents to eliminate; if no restrictions are known set exponents:={}
    // Note that the ouput is {} iff q does NOT bound any exponent.
    // Moreover, if the ouput is different from {}, then it contains q by construction, but this is harmless for the elimination procedure.
    B:=1;
    factQ:=Factorisation(q*OK);
    for x,y in [0..q-1] do
        if ([x,y] ne [0,0]) and (x le y) then
            Bxy:=0; 
            C:=FreyF(x,y);
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
    Sf:={p : p in Set(PrimeDivisors(B)) | p notin {2,3,7,13}};
    if exponents ne {} then
      return (Sf meet exponents) join {q};
    else
      return Sf join {q};
    end if;
end function;


// This function is used to eliminate each form in a given space using specified auxiliary primes and NORMS of differences of traces
function BoundF(forms,AuxiliaryPrimes);
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
                Sq:=BoundFormF(q,f,Sf);
                if Sq ne {} then // Here f can be discarded for large enough p
                    if bool eq 0 then
                        print "This form can be eliminated for large enough p !";
                        Sf:=Sq;
                        bool:=1;
                    end if;
                    Sf:=Sf meet Sq;
                    Sf:={p : p in Sf | p notin {2,3,7,13}};
                end if;
                if Sf ne {} then
                    print "Prime exponent(s) remaining to eliminate =",Sf;
                end if;
            end if;
        end for;
		if bool eq 0 then
		    print "Form no",i," not eliminated for large enough p";
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


// This function is used to eliminate the exponent p for f using a fixed prime q and DIFFERENCES of traces
// Note that the ouput is {} iff q does NOT eliminate p.
// Moreover, if the ouput is different from {}, then it contains q by construction.
// This implies that this function should NOT be used directly. Instead use the next function RefinedBoundF.
function RefinedBoundFormF(q,f,p);
	// q is the auxiliary prime
	// f is the form to eliminate
	// p is the exponent to eliminate
	F:=CoefficientField(f);
	if F eq Rationals() then
		OF:=1;
	else
		OF:=Integers(F);
	end if;
	factQ:=Factorisation(q*OK);
	B:=1*OF;
	for x,y in [0..q-1] do
		if (x le y) and ([x,y] ne [0,0]) then
			Bxy:=0*OF;
			C:=FreyF(x,y);
			for i in [1..#factQ] do
				Q:=factQ[i,1];
				if LocalInformation(C,Q)[3] eq 0 then
					// Here C has good reduction at Q
					diffQ:=TraceOfFrobenius(C,Q) - HeckeEigenvalue(f,Q);
				else
					// Here C has bad multiplicative reduction at Q
					diffQ:=(Norm(Q)+1)^2 - HeckeEigenvalue(f,Q)^2;
				end if;
				if F eq Rationals() then
					Bxy:=Gcd(Integers()!Bxy,Integers()!diffQ);			
				else
					Bxy:=Gcd(Bxy,diffQ*OF);
				end if;
			end for;
			if Bxy eq 0*OF then
            			return {}; // Here p is unbounded
            		else
              			B:=B*Bxy;
            		end if;
		end if;
	end for;
	assert B ne 0*OF;
	Sf:=Set([I[1] : I in Factorisation(q*Gcd(B,p*OF))]);
    return Sf;
end function;

// For a given form f, this function returns the possible ideals P|p (possibly for restricted values of p) such that 
// there is a possible congruence between f and the twisted curve using refined elimination with a given set of auxiliary primes

// This function is used to eliminate the exponent p for f using in turn a bunch of auxiliary primes and DIFFERENCES of traces
// Even in the case where only a single auxilliary prime is necessary, this function should be used instead of the previous one.

function RefinedBoundF(f,AuxiliaryPrimes,p);
	// f is the form to eliminate
	// AuxiliaryPrimes is a set of auxiliary
	// p is the exponent to eliminate
	print "Performing refined elimination to eliminate exponent",p,"with set of auxiliary primes",AuxiliaryPrimes;
	print "";
	F:=CoefficientField(f);
	if F eq Rationals() then
		OF:=1;
	else
		OF:=Integers(F);
	end if;
	Sf:=Set([I[1] :  I in Factorisation(p*OF)]);
  	for q in AuxiliaryPrimes do
		if Sf ne {} then
			print "Dealing with q =",q;
			Sq:=RefinedBoundFormF(q,f,p);
            if Sq ne {} then
                Sf:=Sf meet Sq;
            end if;
            if Sf ne {} then
                print "Prime ideal(s) remaining to eliminate =",Sf;
            end if;
        end if;
	end for;
	if Sf eq {} then
        print "The form is eliminated";
    else
        print "The form with coefficient field :", CoefficientField(f) ;
        print "is not eliminated for",#Sf,"prime ideal(s) above :", p;
    end if;
	print "*************************************************************";
	return "";
end function;





print "";
print "We run standard elimination for the 15 newforms at level Q2*Q3*Q13 using auxilliary primes 5,7,11,17 and 31:";



BoundF(forms111,[5,7,11,17,31]);


print "";
print "The primes p = 5 survives for forms no 6, 9, 10, 13.";
print "The primes p = 11 survives for forms no 6, 9, 10, 13 and 11, 15 respectively.";
print "";
print "We use refined elimination to deal with these remaining forms/exponents";


print "";
print "Dealing with form no 6 and exponent 5:";
RefinedBoundF(forms111[6],[47],5);
// Alternatively:
// RefinedBoundF(forms111[6],[11,31],5);




print "";
print "Dealing with form no 9 and exponent 5:";
RefinedBoundF(forms111[9],[31],5);

print "";
print "Dealing with form no 10 and exponent 5:";
RefinedBoundF(forms111[10],[47],5);

print "";
print "We have discarded exponent 5 for all forms, except for form no 13.";
print "";

print "";
print "Dealing with form no 11 and exponent 11:";
RefinedBoundF(forms111[11],[5],11);

print "";
print "Dealing with form no 15 and exponent 11:";
RefinedBoundF(forms111[15],[5],11);

print "";
print "We have discarded exponent 11 for all forms.";
print "";

print "";
print "This proves Theorem 2, except for p = 5 (form no 13 isn't eliminated for p = 5).";
print "";



// This will finish the case p = 5 for form 13



print "";
print "The unique pair (f,p) which is not eliminated is the pair (f,p) = (form no 13, 5).";
print "We deal with this case.";
print "";





print "Using refined elimination with q = 11, 31 we check there is a unique prime ideal P above 5 which survives.";
f13:=forms111[13]; 
Of13:=Integers(BaseField(f13));

S:={I[1] : I in Factorisation(5*Of13)}; 

print "The size of the residue fields of the primes in Qf13 above 5 are:";
[#ResidueClassField(I) : I in S];


AuxiliaryPrimes:=[11, 31];
for q in AuxiliaryPrimes do 
	S:= S meet RefinedBoundFormF(q,f13,5); 
end for;
assert #S eq 1;
P:=Rep(S);
print "The unique prime P above 5 which we cannot eliminate is:", P;
assert #ResidueClassField(P) eq 5;



print "We check that level raising cannot occur at q=19, so that we can only care about the mod P congruence in the case of good reduction.";
Q:=Factorization(19*OK)[1,1];
assert Valuation(HeckeEigenvalue(f13,Q)^2 - (Norm(Q) + 1)^2,P) eq 0;
print "Done!";


q:=19;
badPairs19:=[];
factQ:=Factorisation(q*OK);
for x,y in [0..q-1] do
    if (x le y) and (x + y) mod q ne 0 then
        Bxy:=0*Of13;
        C:=FreyF(x,y);
        for i in [1..#factQ] do
            Q:=factQ[i,1];
            assert LocalInformation(C,Q)[3] eq 0; 
            diffQ:=TraceOfFrobenius(C,Q) - HeckeEigenvalue(f13,Q);
            Bxy:=Gcd(Bxy,diffQ*Of13);		 
        end for;
        assert Bxy ne 0*Of13;
        if Valuation(Bxy,P) ne 0 then 
            assert Valuation(Bxy,P) gt 0;
            Append(~badPairs19,[x,y]); 
        end if;
    end if;
end for;
print "";
print "The pairs (a,b) mod 19 that make refined elimination fail mod P are ";
badPairs19;



// We now go back to the Frey curve E over Q(\sqrt{13}) 
// From the last sentence of Proposition 9 we know that we have a mod 5 congruence between E and E(1,0), E(1,1) or E(1,-1)
// From the second and third paragraphs in the proof of Theorem 7 we conclude that E is congruent mod 5 to E(1,-1) 
// By taking traces of Frobenius at q=7 using the fact that we know that (a,b) belongs to badPairs7 we obtain a contradiction


print "We now compute the differences between the traces of Frobenius at 19 of E = E(a,b) and E(1,-1) for (a,b) in the list.";

{(TraceOfFrobenius(FreyE(pr[1],pr[2]),19*OF) - TraceOfFrobenius(FreyE(1,-1),19*OF)) : pr in badPairs19};
assert 0 notin $1; 
print "We observe that 0 is not in the previous list, completing the proof.";


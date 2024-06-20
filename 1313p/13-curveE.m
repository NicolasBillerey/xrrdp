/*

We define the elliptic curve E/Q(sqrt(13)) (as in p. 8666 in the published version).

*/

load "13-fieldQsqrt13.m";

L<z13>:=CyclotomicField(13); // Here z13 denotes a primitive 13-th root of unity
RL:=PolynomialRing(L,2);
FL<x,y>:=FieldOfFractions(RL);



f1:= x^2 + (z13 + 1/z13)*x*y + y^2;
f2:= x^2 + (z13^3 + 1/z13^3)*x*y + y^2;
f3:= x^2 + (z13^4 + 1/z13^4)*x*y + y^2;

alpha:= z13^4 + 1/z13^4 - z13^3 - 1/z13^3;
beta:= z13 + 1/z13 - z13^4 - 1/z13^4;
gamma:= z13^3 + 1/z13^3 - z13 - 1/z13;

A:=alpha*f1;
B:=beta*f2;
C:=gamma*f3;

SL<X>:=PolynomialRing(FL);

a4:=3^3*(A*B + A*C + B*C);
a6:=-3^3*(2*A^3 + 3*A^2*B - 3*A*B^2 - 2*B^3);

E:=EllipticCurve(X^3 + a4*X + a6); // The elliptic curve E (but defined over L = Q(zeta13))
AI:=aInvariants(E); // Coefficients of E


/*
// The field F = Q(sqrt(13)) (i.e., the unique quadratic subfield in Q(zeta_13))
F<w>:=QuadraticField(13);
OF:=Integers(F);
*/

RF<x1,y1>:=PolynomialRing(F,2);
_,gm:=IsSubfield(F,L);

AIn:=RL!AI[4]; // a4 coefficient of E
NM:=[Evaluate(c,[x1,y1]) : c in Monomials(AIn)];
NC:=[F!(gm^(-1))(c) : c in Coefficients(AIn)];
NAI4:=RF!(&+[NC[i]*NM[i] : i in [1..#NM]]);

AIn:=RL!AI[5]; // a6 coefficient of E
NM:=[Evaluate(c,[x1,y1]) : c in Monomials(AIn)];
NC:=[F!(gm^(-1))(c) : c in Coefficients(AIn)];
NAI6:=RF!(&+[NC[i]*NM[i] : i in [1..#NM]]);

function FreyE(a,b)
    E:=EllipticCurve([Evaluate(NAI4,[a,b]),Evaluate(NAI6,[a,b])]);
    return E;
end function;




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
                if Sf ne {} then
                    print "Prime exponents remaining to eliminate =",Sf;
                end if;
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
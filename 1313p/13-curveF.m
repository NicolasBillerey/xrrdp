/*

We define the elliptic curve F/Q(sqrt(13)) (as in p. 8669 in the published version).
Here K is the cubic subfield of Q(zeta_13).

*/

L<z13>:=CyclotomicField(13); // Here z13 denotes a primitive 13-th root of unity
RL:=PolynomialRing(L,2);
FL<x,y>:=FieldOfFractions(RL);



f1:= x^2 + (z13 + 1/z13)*x*y + y^2;
f2:= x^2 + (z13^8 + 1/z13^8)*x*y + y^2;

alpha:= z13^8 + 1/z13^8 - z13 - 1/z13;
beta:= 2 - z13^8 - 1/z13^8;
gamma:= z13 + 1/z13 - 2;

A:=alpha*(x+y)^2;
B:=beta*f1;
C:=gamma*f2;

SL<X>:=PolynomialRing(FL);

a4:=3^3*13^2*(A*B + A*C + B*C);
a6:=-3^3*13^3*(2*A^3 + 3*A^2*B - 3*A*B^2 - 2*B^3);

F:=EllipticCurve(X^3 + a4*X + a6); // The elliptic curve F (but defined over L = Q(zeta13))
AI:=aInvariants(F); // Coefficients of F


// The field K (i.e., the unique cubic subfield in Q(zeta_13))
S<X>:=PolynomialRing(Rationals());
K<z>:=NumberField(X^3 + X^2 - 4*X + 1);
OK:=Integers(K);


RK<x1,y1>:=PolynomialRing(K,2);
_,gm:=IsSubfield(K,L);

AIn:=RL!AI[4]; // a4 coefficient of F
NM:=[Evaluate(c,[x1,y1]) : c in Monomials(AIn)];
NC:=[K!(gm^(-1))(c) : c in Coefficients(AIn)];
NAI4:=RK!(&+[NC[i]*NM[i] : i in [1..#NM]]);

AIn:=RL!AI[5]; // a6 coefficient of F
NM:=[Evaluate(c,[x1,y1]) : c in Monomials(AIn)];
NC:=[K!(gm^(-1))(c) : c in Coefficients(AIn)];
NAI6:=RK!(&+[NC[i]*NM[i] : i in [1..#NM]]);

function FreyF(a,b)
    F:=EllipticCurve([Evaluate(NAI4,[a,b]),Evaluate(NAI6,[a,b])]);
    return F;
end function;

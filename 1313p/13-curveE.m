/*

We define the elliptic curve E/Q(sqrt(13)) (as in p. 8666 in the published version).

*/


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


// The field F = Q(sqrt(13)) (i.e., the unique quadratic subfield in Q(zeta_13))
F<w>:=QuadraticField(13);
OF:=Integers(F);


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

I2:=Factorization(2*OF)[1][1];


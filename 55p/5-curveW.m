/*

We define the elliptic curve W/Q and check the computational assertions in Proposition 1.

*/



function FreyW(a,b,d);
    return EllipticCurve([0,-5*(a^2+b^2)*d,0,5*(a^4 - a^3*b + a^2*b^2 - a*b^3 + b^4)*d^2,0]);
end function;

print "";
print "We check that the trace of Frobenius at 3 of the twist by -1 of the elliptic curve W = E(a,b) is -2, 1 or 2:"; 
for a,b in [0..2] do
    if (a le b) and ([a mod 3, b mod 3] ne [0,0]) then
        print TraceOfFrobenius(FreyW(a,b,-1),3);
    end if;
end for;


print "";
print "We check that the trace of Frobenius at 3 of the twist by -1 of the elliptic curve W0 is -1:"; 
W0twisted:=QuadraticTwist(EllipticCurve([0,1,0,592,-16812]),-1);
print "a3(W0^{(-1)}) =",TraceOfFrobenius(W0twisted,3);

print "";
print "This proves Proposition 1.";
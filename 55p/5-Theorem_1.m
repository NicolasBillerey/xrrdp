/*

Here we prove Theorem 1. It uses W/Q (Proposition 1) and E/Q(sqrt(5)) (Part (1) of Proposition 4 and Proposition 5).

*/

load "5-overQsqrt5.m";



print "We apply standard elimination with auxiliary primes 3, 7, 11, 17, and 29 in the case 5 does not divide a + b and p >= 7 using E:";
Bound(Nfs,FreyE,[3,7,11,17,29],7);
print "Done!";
print "There are 6 forms not eliminated.";


print "";
print "We check that these remaining forms correspond to the twists of Frey curves attached to trivial solutions as stated.";

assert SameTraces(Nfs[24],FreyE(1,0,1));
assert SameTraces(Nfs[18],FreyE(1,0,-1));
assert SameTraces(Nfs[17],FreyE(1,0,2));
assert SameTraces(Nfs[23],FreyE(1,0,-2));
assert SameTraces(Nfs[21],FreyE(1,1,1));
assert SameTraces(Nfs[19],FreyE(1,1,2));


print "Done!";
print "This proves part (1) of Proposition 4.";

print "";
print "We finally prove Proposition 5: the traces of Frobenius at 3 is 4 for each elliptic curve E(x,y)^{(d)} that appears in Proposition 4.";
assert {TraceOfFrobenius(curve,I3) : curve in [FreyE(1,0,1), FreyE(1,0, -1), FreyE(1,0,2), FreyE(1,0,-2), FreyE(1,1,1), FreyE(1,1,2)]} eq {4};

print "Done!";
print "Together with Proposition 1 (using W/Q), this proves Theorem 1.";



/*

We check the computational assertions in the proof of Theorem 3 (see Section 6).

*/

Levels:=[50,200,400];


for N in Levels do
    forms:=Newforms(CuspidalSubspace(ModularForms(Gamma0(N),2)));
    print "***********";
    print "There are",#forms,"forms of level",N;
    for i in [1..#forms] do
        f:=forms[i][1];  
        E:=EllipticCurve(f);
        T:=#TorsionSubgroup(E);
        if T mod 2 eq 0 then
            print "The elliptic curve corresponding to the form no",i,"has a rational 2-torsion point.";
            print "";
        end if;
        if T mod 2 eq 1 then
            print "The elliptic curve corresponding to the form no",i,"has no rational 2-torsion point.";
            print "We have a3 =",Coefficient(f,3),"and a7 =",Coefficient(f,7);
            print "";
        end if;
    end for;
end for;

print "This proves Theorem 3.";

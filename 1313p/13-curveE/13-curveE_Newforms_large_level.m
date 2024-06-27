/*

In this file, we compute the newforms of level 2^4*w^2 over Q(sqrt(13)). 
We store their coefficients at the prime ideals in the file 13-curveE_PrimeIdeals.out in the file

13-curveE_Newforms_large_level.out

This is then used to speed up the elimination procedure at this level in the proof of Proposition 9.

*/


F<w>:=QuadraticField(13);
OF:=Integers(F);
I2:=2*OF;
I13:=Factorisation(13*OF)[1,1];
I3a:=Factorisation(3*OF)[1,1];
I3b:=Factorisation(3*OF)[2,1];
R<e>:=PolynomialRing(Rationals());






PrimeIdeals:=[];
for pi in PrimeGenerators do
    Append(~PrimeIdeals,(OF!pi)*OF);
end for;

load "13-curveE_PrimeIdeals.out";

procedure PrintNewforms(filename,label,forms);

    PrintFile(filename, label cat ":=[":Overwrite:=true);

    for i:=1 to #forms do
        f:=forms[i];

        EI:=[];
        for Q in PrimeIdeals do
            Append(~EI,HeckeEigenvalue(f,Q));
        end for;

        KF:=Parent(EI[1]);
        FI:=<R!DefiningPolynomial(KF),EI>;

        PrintFile(filename,FI);
        if i lt #forms then
            PrintFile(filename,",");
        end if;
    end for;

    PrintFile(filename, "];");

end procedure;



N2:=I2^4*I13^2;
print "Computing newforms of level 2^4*w^2. Space of dimension", Dimension(NewSubspace(HilbertCuspForms(F,N2)));
time NewformsN2:=Eigenforms(NewSubspace(HilbertCuspForms(F,N2)));
print "Done!";
PrintNewforms("13-curveE_Newforms_large_level.out","CoefficientsN2",NewformsN2);
# A multi-Frey approach to Fermat equations of signature (r,r,p)

!!!CAREFUL: WORK IN PROGRESS!!!

Electronic resources for the paper `A multi-Frey approach to Fermat equations of signature (r,r,p)' (Trans. Amer. Math. Soc. 371 (2019), no. 12, 8651--8677) and its <a href="https://arxiv.org/abs/1703.06530">arXiv post-publication revised version</a> by Nicolas Billerey, Imin Chen, Luis Dieulefait, and Nuno Freitas.

Remark: The programs were run on a 2.35/3.35 Ghz 32 core AMD EPYC 7452 machine with 512 Gb from Laboratoire de Mathématiques Blaise Pascal in Université Clermont Auvergne using Magma V2.28-9.

Last modifications: Jun 18, 2024

********************************
Part I - signature (5,5,p)
********************************

5-curveW.m : definition of W/Q and proof of Proposition 1.

5-overQsqrt5.m : useful definitions (such as E and F over Q(sqrt5)) and space computation (for Proposition 4)

5-Theorem_1.m : proof of Theorem 1 (and hence of Part (1) of Proposition 4)

5-Theorems_5_and_6.m : proof of Theorems 5 and 6 (and hence of Part (2) of Proposition 4)

5-SecondCase.m : proof of Theorem 3 (see Section 6)


********************************
Part II - signature (13,13,p)
********************************

* Files related to the curve E:

13-curveE.m : definition the elliptic curve E/Q(sqrt(13)) (this definition has been updated to match the definition of the paper; the previous code was using a different model of this curve)

13-curveE_elimination_large_level_first_approach.m : proof of Proposition 9 in the "large" level case using precomputation of newforms at that level 2^4*w^2 (from files 13-curveE_Newforms_large_level.out and 13-curveE_PrimeIdeals.out)

13-curveE_elimination_large_level_second_approach.m : proof of Proposition 9 in the "large" level case using classical elimination (STILL COMPUTING!)

13-curveE_elimination_small_level.m : proof of Proposition 9 in the "small" level case using classical elimination

13-curveE_inertia.m : proof of (part (B) of) Theorem 7

13-curveE_irreducibility.m : proof of Proposition 8

13-curveE_Newforms_large_level.m: compute the list of coefficients of forms of level 2^4*w^2 and store them in 13-curveE_Newforms_large_level.out (DOES NOT WORK WELL)

13-curveE_PrimesIdeals.m: compute a list of prime ideals of Q(sqrt(13)) of norm up to 200 and store them in the file 13-curveE_PrimeIdeals.out.

* Files related to the curve F:

13-curveF.m : definition the elliptic curve F/Q(sqrt(13))

13-curveF_elimination.m : proof of Theorem 2

13-curveF_irreducibility.m : proof of Lemma 10



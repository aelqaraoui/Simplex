# Simplex

-0.5 -3 -1 -4 0 0 0 0 0
1 1 1 1 1 0 0 0 40
2 1 -1 -1 0 1 0 0 10
0 -1 0 1 0 0 0 1 10

n = 9, m = 4, a = 0, ukn = 4












Maximize z = - x1 - x2 + 0s1 + 0s2 - A1 - A2
subject to the constraints :
2x1 + x2 – s1 + A1 = 4,
x1 + 7x2 – s2 + A2 = 7,
and x1 , x2, s1, s2 , A1, A2 ≥ 0. 

converted to give matrix input
also input number of artificial vars 

Maximize z = - x1 - x2 + 0s1 + 0s2 - A1 - A2
subject to the constraints :
2x1 + x2 – s1 +0s2 + A1 + 0A2 = 4,
x1 + 7x2 +0s1 – s2 + 0A1 + A2 = 7,
and x1 , x2, s1, s2 , A1, A2 ≥ 0. 

assuming all constraints are of form Ax<=b , where b can be negative.(so no equality constraint)

corresponding matrix in std input is
//note z all shifted left 
1 1 0 0 1 1 0      		
2 1 -1 0 1 0 4     
1 7 0 -1 0 1 7

feasible

-------------------------------------------------------------------------------------------------------
Maximize z = x1 - 2x2 – 3x3 - A1 - A2
subject to the constraints :
-2x1 + x2 + 3x3 + A1 = 2,
2x1 + 3x2 + 4x3 + A2 = 1,
and x1 , x2, A1, A2 ≥ 0

numartific = 2
slacks = 0

Maximize z - x1 + 2x2 + 3x3 + A1 + A2 = 0
subject to the constraints :
-2x1 + x2 + 3x3 + A1 + 0A2 = 2,
2x1 + 3x2 + 4x3 + 0A1 + A2 = 1,
and x1 , x2, A1, A2 ≥ 0


-1 2 3 1 1 0
-2 1 3 1 0 2
2 3 4 0 1 1

Not feasible
--------------------------------------------------------------------------------------------------------------
Maximize Z = 2x1 + 4x2,
subject to 2x1 + 3x2 ≤ 48,
 x1 + 3x2 ≤ 42,
 x1 + x2 ≤ 21,
 and x1, x2 ≥ 0 

Maximize Z = 2x1 + 4x2 + 0s1 + 0s2 + 0s3,
subject to 2x1 + 3x2 + 1s1 ≤ 48,
 x1 + 3x2 +1s2 ≤ 42,
 x1 + x2 + 1s3 ≤ 21,
 and x1, x2 ≥ 0

-2 -4 0 0 0 0
2 3 1 0 0 48
1 3 0 1 0 42
1 1 0 0 1 21

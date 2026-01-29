We want to analyse the evolution of of sales of an enterprise .
## 1-The problem can be write like a linear model  :
 min ||A.theta-b||(minimization problem)

 with A[] = {
    1.0, 10.0, 3.0,1.2,
    1.0,12.0,4.0,1.1,
    1.0,8.0,6.0,0.9,
    1.0,15.0,5.0,1.3,
    1.0,9.0,6.0,1.0
}
and b[] = {25.0, 28.0,20.0, 32.0, 22.0}

## 2-The solution we obtain with the lapack function LAPACKE_dgels is theta = [11.574194, 1.677419, -0.612903, -1.225806]
so a0=11.574194;
   a1=1.677419;
   a2=-0.612903;
   a3=-1.225806;

## 3- Interpretation
a0 is the bias of the model
a1 ,a2 and a3 are the weights those minimize the gap between the model approximation and the real sale prices.

## 4-Conclusion 


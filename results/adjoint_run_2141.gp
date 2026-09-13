
default(parisize, "1G"); default(realprecision, 30);
read("/home/claude/kit_even_icosahedral/scripts/../results/adjoint_coefficients_2141.gp");
Ncond = 2141^2;
{
for(i=1,2,
  eps = if(i==1, 1, -1);
  L = lfuncreate([Ad_an, 0, [0,0,0], 1, Ncond, eps]);
  print("eps = ", eps, ":  lfuncheckfeq (log2 of FE residual) = ", lfuncheckfeq(L));
);
L = lfuncreate([Ad_an, 0, [0,0,0], 1, Ncond, 1]);
L1 = lfun(L, 1); print("L(1, Ad) = ", L1);
print("L(1/2, Ad) = ", lfun(L, 1/2));
print("L(2, Ad) = ", lfun(L, 2));
print("predicted RS constant L(1,Ad)/zeta(2)*2141/2142 = ", L1/zeta(2)*2141/2142);
}

default(parisize,"512M");
N=2141; mu = N+1; e2 = 1+kronecker(-1,N); e3 = 1+kronecker(-3,N); c=2;
print("2141 prime? ", isprime(N), " index=",mu," e2=",e2," e3=",e3," genus=", 1 + mu/12 - e2/4 - e3/3 - c/2);
P = mspolygon(N); E = P[3];
f = fileopen("generators_2141.txt", "w");
for(i=1,#E, M=E[i]; if(type(M)=="t_MAT", filewrite(f, Str(M[1,1]," ",M[1,2]," ",M[2,1]," ",M[2,2]))));
fileclose(f); print("wrote ", #E);

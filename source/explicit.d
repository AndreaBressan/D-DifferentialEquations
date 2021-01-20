/**
* Module for templated explicit steps for ODE solvers.
* The ODE is assumed to be in canonical form:
*
* y'(t)=f(t,y(t))
*
* Given the t, y(t) and a time increment dt a stepping strategy produces
* an approximation of y(t+dt).
*
* Each step is a templated function that takes 3 or more parameters:
* * the function f
* * the starting time t
* * the time increment dt
* and it returns y(t+dt)
*
* The function are templated so that arbitrary types can be used for
* the value f(t) allowing the solution of scalar-, vector- or tensor-valued
* ODEs and different implementations of scalars.
*
* All the method are special cases of Runge-Kutta and they share a generic
* implementation that is based on the Butcher table of the coefficients of
* the method.
*/
module explicit.d;

import std.traits: Unqual, hasMember;
import std.typecons;

immutable simpleStepNames = ["euler", "rk3", "rk4","heun"];
immutable adaptiveStepNames = ["dopri"];


/// This loop generates the functions for the different methods
/// that return the computed step
static foreach(name; simpleStepNames)
{
mixin("valT "~name~"Step (funT, timeT, valT)( funT f, timeT t, timeT dt, valT y){return rkStep!("~name~"Table,funT,timeT,valT)(f,t,dt,y);} " );
}

/// This loop generates the functions for the different methods
/// that return to the computed with two different methods and
/// allow for error estimation
static foreach(name; adaptiveStepNames)
{
mixin("Tuple!(Unqual!valT,Unqual!valT) "~name~"Step (funT, timeT, valT)( funT f, timeT t, timeT dt, valT y){return rkStep!("~name~"Table,funT,timeT,valT)(f,t,dt,y);} " );
}

/// Explicit Runge-Kutta methods are a multi-step method where intermediate
/// results are combined according to a predetermined Butcher table.
auto rkStep(alias table,funT, timeT, valT)( funT f, timeT t, timeT dt, valT y)
/// TODO add template constraints
{
/// util to be factored out, modeled over std.numeric.dotProduct
/// that cannot be used as it assumes that a and b are both scalars
Unqual!valT weightedSum(const Unqual!valT[] v, const Unqual!timeT[] s) pure
{
// avoid the problem of having a zero inizializer
Unqual!valT result = s[0]*v[0]; 
for (auto i=1;i<s.length;++i)
{
result += s[i] * v[i];
}
return result;
}

alias a=table!(timeT).a;
alias b=table!(timeT).b;
alias c=table!(timeT).c;

immutable len=a.length;
Unqual!valT[len] k;
k[0]=f(t,y);
foreach ( i; 0.. len)
{
k[i]=f(t+c[i]*dt, weightedSum(k[0..i],a[i][0..i]));
}
static if (hasMember!(table!(timeT),"b2"))
{
    alias b2=table!(timeT).b;
    return tuple(weightedSum(k,b),weightedSum(k,b2));
}
else
    return weightedSum(k,b);
}



struct rk3Table(T)
{
static immutable a=[
[T(0)],
[T(0), T(1)/T(2)],
[T(0), T(0), T(1)/T(2)]];
static immutable b=[T(1)/T(6),T(1)/T(3),T(1)/T(3),T(1)/T(6)];
static immutable c=[T(0),T(1)/T(2),T(1)/T(2),T(1)];
}

struct rk4Table(T)
{
static immutable a=[
[T(1)/T(2)],
[T(0), T(1)/T(2)],
[T(0), T(0), T(1)/T(2)]];
static immutable b=[T(1)/T(6),T(1)/T(3),T(1)/T(3),T(1)/T(6)];
static immutable c=[T(0),T(1)/T(2),T(1)/T(2),T(1)];
}


struct eulerTable(T)
{
static immutable a=[[T(0)]];
static immutable b=[T(1)];
static immutable c=[T(0)];
}

struct heunTable(T)
{
static immutable a=[
[T(0)],
[T(1), T(0)]];
static immutable b=[T(1)/T(2), T(1)/T(2)];
static immutable c=[T(0), T(1)];
}

struct dopriTable(T)
{
static immutable a=[
[T(1)/T(5)],
[T(3/40), T(9)/T(40)],
[T(44)/T(45), T(-56)/T(15), T(32)/T(9)],
[T(19372)/T(6561), T(-25360)/T(2187), T(64448)/T(6561), T(-212)/T(729)],
[T(9017)/T(3168), T(-355)/T(33), T(46732)/T(5247), T(49)/T(176), T(-5103)/(18656)],
[T(35)/T(384), T(0), T(500)/T(1113), T(125)/T(192), T(-2187)/T(6784), T(11)/T(84)]
];
static immutable b=[T(5179)/T(57600), T(0), T(7571)/T(16695), T(393)/T(640), T(-92097)/T(339200), T(187)/T(2100),T(1)/T(40)];
static immutable c=[T(0),T(1)/T(5),T(3)/T(10),T(4)/T(5), T(8)/T(9), T(1), T(1)];
static immutable b2=[T(0), T(0), T(0), T(0), T(0), T(0), T(1)];
}



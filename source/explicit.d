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
module explicit.d

import std.traits:Unqual;

/// This loop generates the functions for the different methods
/// from a list of names
static foreach(name; ["euler", "rk3", "rk4","heun"])
{
  mixin("valT "~name~"Step (funT, timeT, valT)( funT f, timeT t, timeT dt, valT y){return rkStep!("~name~"Table,funT,timeT,valT)(f,t,dt,y);} " );
}

/// Explicit Runge-Kutta methods are a multi-step method where intermediate
/// results are combined according to a predetermined Butcher table.
Unqual!valT rkStep(alias table,funT, timeT, valT)( funT f, timeT t, timeT dt, valT y)
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

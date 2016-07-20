
module PolyaGamma

using StatsBase, StatsFuns

export rpolyagamma

function rpolyagamma(z::Float64,t::Float64=0.64)
  accept = false;
  z = abs(z)/2;
  K = pi^2/8 + z^2/2;
  p = pi/(2K) * exp(-K*t);
  q = 2*exp(-z)*pinvgauss(t,1/z);
  x::Float64 = 0.0;

  while !accept
    #generate proposal
    if rand() < p/(p+q)
      x = t + randexp()/K;
    else
      x = rtinvgauss(1/z,t);
    end

      #iterate sum until accept/reject
    reject = false;
    n = 0;
    s = asterm(x,n);
    y = s*rand();
    while !accept & !reject
      n += 1;
      if isodd(n)
        s -= asterm(x,n);
        accept = y < s;
      else
        s += asterm(x,n);
        reject = y > s;
      end
    end
  end

  return x/4.0
end

function rpolyagamma(z::Float64,b::Int64,t::Float64=0.64)
  x = 0.0;
  for i in 1:b x += rpolyagamma(z,t); end
  return x
end


function pinvgauss(x::Float64,μ::Float64,λ::Float64=1.0)
  z = sqrt(λ/x) * (x/μ - 1);
  zinv = -sqrt(λ/x) * (x/μ + 1);
  return normcdf(z) + exp(2λ/μ) * normcdf(zinv);
end

function asterm(x,n::Int64,t::Float64=0.64)
  pin = pi*(n+0.5);
  if x <= t
    a = pin * (2/(pi*x))^1.5 * exp(-2*(n+0.5)^2/x);
  else
    a = pin * exp(-pin^2*x/2);
  end
  return a
end

function rtinvgauss(μ::Float64,t::Float64 = 0.64)
  if μ > t

    a = 0.0;
    while rand() > a
        e1 = randexp(); e2 = randexp();
      while e1^2 > (2*e2/t)
        e1 = randexp(); e2 = randexp();
      end
      x = t/(1+t*e1)^2;
      a = exp(-x/(2*μ^2));

    end

  elseif μ <= t

    x = t+1.0;
    while x > t
      y = randn()^2;
      x = μ + 0.5*y*μ^2 - 0.5*μ*sqrt(4*μ*y + μ^2*y^2);
      if rand() > μ/(μ+x)
        x = μ^2/x;
      end
    end

  end

  return x
end

end

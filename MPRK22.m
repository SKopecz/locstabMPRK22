function [t, Y] = MPRK22(y0, dt0, tint, a21, a, b)

assert(a21 ~= 0)

y0 = y0 + realmin;

b2 = 1/(2*a21);
b1 = 1 - b2;

p = 1/a21;
q = 1 - p;

I = length(y0);

tstart = tint(1);
tend = tint(2);

t = tstart;
dt = dt0;

P = zeros(2);
D = P;
P2 = P;
D2 = P;

Y = y0;
y = y0;
while t(end) < tend
  if t(end) + dt > tend
    dt = tend - t(end);
  end

  %[P, D] = PD(t(end),y);
  by2 = b*y(2);
  ay1 = a*y(1);
  P(1,2) = by2;
  P(2,1) = ay1;
  D(1,2) = ay1;
  D(2,1) = by2;

  M2 = zeros(I);
  for i = 1:I
    if a21 > 0
      M2(i,i) = 1 + a21*dt*sum(D(i,:))/y(i);
    else
      M2(i,i) = 1 - a21*dt*sum(P(i,:))/y(i);
    end
    for j = 1:I
      if j ~= i
        if a21 > 0
          M2(i,j) = -a21*dt*P(i,j)/y(j);
        else
          M2(i,j) =  a21*dt*D(i,j)/y(j);
        end
      end
    end
  end
  y2 = M2\y;

  y2 = y2 + realmin;

  %[P2, D2] = PD(t(end)+a21*dt,y2);
  by22 = b*y2(2);
  ay21 = a*y2(1);
  P2(1,2) = by22;
  P2(2,1) = ay21;
  D2(1,2) = ay21;
  D2(2,1) = by22;

  M = zeros(I);
  sigma = ((y2./y).^p).*y;
  for i = 1:I
    if b1 >= 0
      M(i,i) = 1 + dt*b1*sum(D(i,:))/sigma(i);
    else
      M(i,i) = 1 - dt*b1*sum(P(i,:))/sigma(i);
    end
    if b2 >= 0
      M(i,i) = M(i,i) + dt*b2*sum(D2(i,:))/sigma(i);
    else
      M(i,i) = M(i,i) - dt*b2*sum(P2(i,:))/sigma(i);
    end
    for j=1:I
      if j ~= i
        if b1 >= 0
          M(i,j) = -dt*b1*P(i,j)/sigma(j);
        else
          M(i,j) = dt*b1*D(i,j)/sigma(j);
        end
        if b2 >= 0
          M(i,j) = M(i,j) - dt*b2*P2(i,j)/sigma(j);
        else
          M(i,j) = M(i,j) + dt*b2*D2(i,j)/sigma(j);
        end
      end
    end
  end
  y = M\y;

  y = y + realmin;

  Y(:,end+1) = y;
  t(end+1) = t(end) + dt;
end
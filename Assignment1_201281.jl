##Functions Declaration
function GradDesc (xInit, yPrime, precision, eps)
  xOld=xInit;
  xNew=6;
  iterations=0;
  histX=[iterations];
  while abs(xNew - xOld)>precision
    xOld=xNew;
    xNew=xOld- eps*yPrime (xOld);
    push!(histX , iterations+1)
    
  end

  print ("[Solution]: ", xNew);
  xNew;
  print ("    |||     [Steps]: ",length(histX))

end
##################
function NewtMeth(xInit, yPrime,yDPrime, precision)
  xOld=0;
  xNew=3/10;
  iterations=0;
  histX=[iterations];
  
  while abs(xNew - xOld)>precision
    xOld=xNew;
    xNew=xOld- yPrime(xOld) / yDPrime(xOld);
    push!(histX , iterations+1)
    
  end
  print ("[Solution]: ", xNew);
  xNew;
  print ("     |||    [Steps]: ",length(histX))
  
end
#################
xInit=0;
precision=0.00001;
eps=0.00001;


#1st Function [Minimize]
#y(x)=1.2*(x-2)^2 + 3.2
#yPrime(x)=12/5 * (x-2);
#yDPrime(x)= 12/5

#2nd Function [Minimize]
#y(x)=ln (1+x^2)
#yPrime(x)=2*x/ (x^2 +1);
#yDPrime(x)= (2/ (x^2 +1))-(4*x^2 / (x^2 +1 )^2)

#3rd Function [Maximize]
#y(x)=-x^2 + 2*x;
#yPrime(x)=2*x-2
#yDPrime(x)=2
function initpopReal(npop, R_max, R_min)

  initialPopReal1=rand(R_min[1]:R_max[1], npop,1)
  initialPopReal2=rand(R_min[2]:R_max[2], npop,1)

  initialPopReal=hcat(initialPopReal1,initialPopReal2)
  return initialPopReal
end

function evalFitnessPop (pop)

  popFit=zeros(Int64,1,npop)

  for i=1:npop
    fitness=8-(pop[i,1])^2 - (pop[i,2])^2
    if i==1
      popFit=vcat(fitness)
    elseif (i>1)
      popFit=vcat(popFit,fitness)
    end
  end
  return popFit
end
##########################3
function selectionProb(popFit)
  counter=1
  totalFitness=sum(popFit)
  probs=zeros(length(popFit),1)
  while(counter<=length(popFit))
    probs[counter]=popFit[counter]/totalFitness
    counter+=1
    end
  return probs
end

function cum_prob(probs)
    cprob = cumsum (probs)
    return cprob
end

function roulette_wheel(cprob)
    R=rand()
    for index=1:length(cprob)
      if(cprob[index]>R)
        return (index)
      end
    end
end

function RealRoulette_select(cprob,pop)
    parent1=roulette_wheel(cprob)
    parent2=roulette_wheel(cprob)

    two_parents=zeros(2,2)
    two_parents[1,:]=pop[parent1,:]
    two_parents[2,:]=pop[parent2,:]
    return two_parents
end

function  arithmetic_crossover (two_parents, pcross)
  two_children=zeros(2,2)
  crossProbability=rand()
  if (crossProbability<pcross)
    for i=1:2
      w=rand()
      two_children[i,1]=w*two_parents[1,1]+(1-w)*two_parents[2,1]
      two_children[i,2]=w*two_parents[2,1]+(1-w)*two_parents[2,2]
    end
  else
    for i=1:2
      two_children[i,1]=two_parents[i,1]
      two_children[i,2]=two_parents[i,2]
    end
  end

return two_children
end

function  gaussian_mutate(individual,sigma,pmute,R_max, R_min)

  m=zeros(1,length(individual))
  mutationr=rand()

  if (mutationr<=pmute)
    for i=1:2
      mean=0
      u1 = rand()
      u2 = rand()
      r = sqrt( -2.0*log(u1) )
      theta = 2.0*pi*u2
      m = mean + sigma*r*sin(theta)
      individual=individual+m
    end

    for j=1:2
      if individual[1,j]>R_max[1,j]
        individual[1,j]=R_max[1,j]

      elseif individual[1,j]<R_min[1,j]
        individual[1,j]=R_min[1,j]
      end
    end

  end
  mutatedIndividual=individual

  return mutatedIndividual
end

function  runRealGA(npop,ngen,pcross,pmute,sigma,R_max,R_min)

  pop=zeros(Int64,npop,2);
  pop=initpopReal(npop,R_max,R_min);
  best_hist=zeros(ngen,1);
  mutatedChildren=zeros(2,2);
  generation= zeros(npop,2);

  for z=1:ngen
    popFit=evalFitnessPop(pop);
    probs=selectionProb(popFit);
    if (z==1)
      best_hist=vcat(maximum(popFit));
    elseif (z>1)
      best_hist=vcat(best_hist,maximum(popFit));
    end
    comulativeProbs=cum_prob(probs);
    for k=1:50
      twoParents=roulette_select(comulativeProbs,pop);
      twoChildren=arithmetic_crossover(twoParents,pcross);
      for i=1:2
        mutatedChildren[i,:]=gaussian_mutate(twoChildren[i,:],sigma,pmute,R_max,R_min)
      end
      if (k==1)
        generation=vcat(mutatedChildren)
      elseif (k>1)
        generation=vcat(generation,mutatedChildren)
      end
    end
    pop=generation
  end
  finpop= pop
  return finpop, best_hist
end

npop=100;
ngen=200;
pcross=0.7;
pmute = 0.01;
sigma = 0.8;
R_max = [2 2];
R_min = [-2 -2];

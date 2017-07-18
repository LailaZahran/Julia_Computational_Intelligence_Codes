# Binary GA (Onemax function):
function initPopBinary(npop,clen)
  initialPop=rand(0:1, npop,clen)
  return initialPop
end

function oneMaxI(chromosome)
  chromosomeFitness=sum(chromosome)
  return chromosomeFitness
end

function oneMaxP(pop,npop)
  counter=1
  popFit=zeros(npop)
  while(counter<=npop)
    popFit[counter]=oneMaxI(pop[counter,:])
    counter+=1
  end
  return popFit
end

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

function roulette_select(cprob,pop)
    parent1=roulette_wheel(cprob)
    parent2=roulette_wheel(cprob)

    two_parents=zeros(Int64, 2,clen)
    two_parents[1,:]=pop[parent1,:]
    two_parents[2,:]=pop[parent2,:]
    return two_parents
end

function bin_cross(two_parents, pcross)
  crand=rand();
  two_children=zeros(2,clen)
  crossPoint = 3
  if (crand <= pcross)
    two_children[1,:] = [two_parents[1,1:crossPoint] two_parents[2,crossPoint+1:end]]
    two_children[2,:] = [two_parents[2,1:crossPoint] two_parents[1,crossPoint+1:end]]

  else
    two_children[1,:]=two_parents[1,:]
    two_children[2,:]=two_parents[2,:]
  end
    return two_children
end

function bin_mutate(individual,pmute)

  for i=1:clen
    rNum=rand();
    if (rNum<pmute)
      if individual[i]==0
        individual[i]=1
      else
        individual[i]=0
      end
    else
    end
  end
  mutated_ind=zeros(Int64,1,clen)
  mutated_ind=individual
  return mutated_ind
end

function runBinGA(npop, clen, ngen, pcross, pmute)

      pop=zeros(Int64,npop,clen);
      pop=initPopBinary(npop,clen);
      mutatedChildren=zeros(Int64,2,clen);
      generation= zeros(Int64,npop,clen);
      best_hist=zeros(Int64,ngen,1);

      for z=1:ngen
        popFit=oneMaxP(pop,npop);
        probs=selectionProb(popFit);
        if (z==1)
          best_hist=vcat(maximum(popFit));
        elseif (z>1)
          best_hist=vcat(best_hist,maximum(popFit));
        end
        comulativeProbs=cum_prob(probs);
        for k=1:50
          twoParents=roulette_select(comulativeProbs,pop);
          twoChildren=bin_cross(twoParents,pcross)

             for i=1:2
                 mutatedChildren[i,:]=bin_mutate(twoChildren[i,:],pmute)
             end
          if (k==1)
           generation=vcat(mutatedChildren)
         elseif (k>1)
        generation=vcat(generation,mutatedChildren)
      end
       end
       pop=generation;
     end

     finpop= pop
return finpop,int(best_hist)
end

      npop = 100;
      clen = 20;
      ngen = 200;
      pcross = 0.7;
      pmute = 0.001;
      result=runBinGA(npop, clen, ngen, pcross, pmute)


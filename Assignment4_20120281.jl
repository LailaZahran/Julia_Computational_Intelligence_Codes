#Computational Intelligence
#Assignment #4
#Ant Colony Optimization


function readCityCoordinates(fileName)
  cd (dirname(@__FILE__))
  citiesCoordinates = readcsv (fileName, Float64)
  numOfCities=size(citiesCoordinates,1)

  return citiesCoordinates,numOfCities
end

function distCalculate(citiesCoordinates,numOfCities)
  distanceBetweenCities=zeros(numOfCities,numOfCities)

  for i= 1:size(citiesCoordinates,1)
    for j=1:size(citiesCoordinates,1)

      distanceBetweenCities[i,j]=sqrt((citiesCoordinates[j,1]-citiesCoordinates[i,1])^2 + (citiesCoordinates[j,2]-citiesCoordinates[i,2])^2)

    end
  end
  return distanceBetweenCities
end

function etaCalculate(distanceBetweenCities,numOfCities)

  eta=zeros(numOfCities, numOfCities)

  for i=1:numOfCities
     for j=1:numOfCities
      eta[i,j]= 1 / (distanceBetweenCities[i,j])
    end
  end

  return eta
end

function heuristicCalculate(distanceBetweenCities,numOfCities)

  dnew=zeros(size(distanceBetweenCities,1), size(distanceBetweenCities,2))
  for i=1:numOfCities
    for j=1:numOfCities
      dnew[i,j]=distanceBetweenCities[i,j]
      end
  end

  for i=1:numOfCities
    for j=1:numOfCities
     if (i==j)
        dnew[i,j]=10000
      end
    end
  end

  Lnn=0
  indexValues=0
  #Start With City#1
  fin=1
  in=fin
  endpoint=0
  indexValues=hcat(in)

  for i=1:numOfCities
    if (i<numOfCities)
      len, to=findmin(dnew[in,:])
      Lnn=Lnn+len
      indexValues=hcat(to,indexValues)
      dnew[:,in]=1000
      in=to

    elseif (i==numOfCities)
      len=distanceBetweenCities[in,fin]
      to=fin
      Lnn=Lnn+len
      in=to
      endpoint=in
    end
  end
  tau0 = 1/(numOfCities*Lnn)
  #Reverse an array
  visitedCities=indexValues[end:-1:1,end:-1:1]

  return  Lnn, tau0
end

function constructPhermonesMatrix(tau0,numOfCities)
  tau=zeros(numOfCities,numOfCities)
  for i=1:numOfCities
    for j=1:numOfCities
      tau[i,j]=tau0
    end
  end
  return tau
end

function stateTransition(visitedCities, tau, eta, q0, beta)

  q=rand()
  numOfCities=size(tau,1)
  fromCity=visitedCities[end]
  state= tau[fromCity,:].* (eta[fromCity,:].^beta)
  toCity=0

  if q<q0
    max=0
    for i=1:numOfCities
      if issubset([i],visitedCities)==false && i!=fromCity
        if state[1,i]>max
          max=state[1,i]
          toCity=i
        end
      end
    end
    nextCity=toCity
    return nextCity

  else
    #Using roulette wheel
    probs=zeros(numOfCities)
		cprobs=zeros(numOfCities)
    total=0
    #Calculate total probabilities
    for i=1:numOfCities
			if( issubset([i],visitedCities)==false && i!=fromCity )
				total+=state[i]
			end
		end
    #Calculate prob. of each city
    for i=1:numOfCities
			if( issubset([i],visitedCities)==false && i!=fromCity )
				probs[i]=state[i]/total
			end
		end
    cprob=cumsum(probs)
		r=rand()
		for i=1:numOfCities
			if( issubset([i],visitedCities)==false && i!=fromCity )
				if(r < cprob[i])
					toCity=i
				end
			end
		end
    nextCity=toCity
    return nextCity
  end
end

function updateLocalPheromone(fromCity, toCity, arcTau, tau0, rho)

  updatedArcTau = (1-rho)*arcTau + rho*tau0

  return updatedArcTau
end

function updateGlobalPhermone(fromCity, toCity, lenBestTour, arcTau,alpha)

  updatedArcTau= (1-alpha)*arcTau + alpha*(1/ lenBestTour)

  return updatedArcTau
end

function runACO(ngen,beta,rho,alpha,q0,numAnts,citiesCoordinates,numOfCities)
  # ngen: number of generations
  #Initialization
  euclideanDistance= distCalculate(citiesCoordinates,numOfCities)
  eta=etaCalculate(euclideanDistance,numOfCities)
  Lnn,tau0=heuristicCalculate(euclideanDistance,numOfCities)
  tau=constructPhermonesMatrix(tau0,numOfCities)
  lenBestTour=1000
  theBestTour=zeros(1,numOfCities)

  for a=1:ngen
####################
  citieschosen=0
  positions=zeros(numAnts,1)
  for i=1:size(positions,1)
    if (i==1)
      positions [i,:]=rand(1:29)
      citieschosen=hcat(positions [i,:],citieschosen)
    else
      positions [i,:]=rand(1:29)
      while issubset(positions [i,:] ,citieschosen)
        positions [i,:]=rand(1:29)
      end
      citieschosen=hcat(positions [i,:],citieschosen)
    end
  end

  antsTour=zeros(numAnts,numAnts)
  antTourLength=zeros(numAnts,1)

  for i=1:size(antsTour,1)
    antsTour[i,1]=positions[i,1]
  end

  for i=1:numAnts

      for j=1:numAnts-1
        currentCity=antsTour[i,j]
        nextCity=stateTransition(antsTour[i,1:j],tau,eta,q0,beta)
        antTourLength[i]=antTourLength[i] + euclideanDistance[currentCity,nextCity]
        tau[currentCity,nextCity]=updateLocalPheromone(currentCity, nextCity, tau[currentCity,nextCity], tau0, rho)
        antsTour[i,j+1]=nextCity
      end
    tau[antsTour[i,end],antsTour[i,1]]=updateLocalPheromone(antsTour[i,end],antsTour[i,1], tau[antsTour[i,end],antsTour[i,1]], tau0, rho)
    antTourLength[i]=antTourLength[i]+ euclideanDistance[antsTour[i,end],antsTour[i,1]]

  end

  bestTour=size(1,numOfCities)
  genlenBestTour,antNumber=findmin(antTourLength)
  bestTour=copy(antsTour[antNumber,:])

  for i=1:numOfCities-1
    tau[bestTour[i],bestTour[i+1]]= updateGlobalPhermone(bestTour[i],bestTour[i+1], genlenBestTour, tau[bestTour[i],bestTour[i+1]],alpha)
  end

  tau[bestTour[end],bestTour[1]]=updateGlobalPhermone(bestTour[end],bestTour[1], genlenBestTour, tau[bestTour[end],bestTour[1]], alpha)

    if (genlenBestTour<lenBestTour)
      lenBestTour=genlenBestTour
      theBestTour=antsTour[antNumber,:]
    end

  end
  ####################
  return lenBestTour, theBestTour
end

citiesCoordinates,numOfCities = readCityCoordinates("TSPdata.csv")
lenBestTour, theBestTour = runACO(100,2,0.1,0.1,0.9,29,citiesCoordinates,numOfCities)

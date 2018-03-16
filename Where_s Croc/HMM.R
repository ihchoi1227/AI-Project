numWaterpools = 40

getAdjacencyMatrix=function(edges) {
	adj = matrix(0, nrow=numWaterpools, ncol=numWaterpools)
	# mark edges in one direction
	adj[edges] = 1
	# mark edges in other direction (symmetrically)
	adj = adj + t(adj)
	# mark each waterpool as adjacent to itself
	for (i in 1:numWaterpools) {
		adj[i,i] = 1
	}
	return (adj)
}

initialState=function(positions) {
	stateVector = rep(1, numWaterpools)
	stateVector = stateVector / sum(stateVector)
	return (stateVector)
}

# All actions from a state are equiprobable, so the
# transition matrix is just the adjacency with normalized rows.
# Rows are normalized because state vector is a row vector,
# so to get next state vector we multiply previous one by matrix.
transitionMatrix=function(adj) {
	for (i in 1:numWaterpools) {
		adj[i,] = adj[i,] / sum(adj[i,])
	}
	return (adj)
}

distanceMatrix=function(adj) {
	distances = adj
	# set distance between yet disconnected node to infinity
	distances[which(distances == 0)] = Inf

	temp = adj

	step = 2
	while (Inf %in% distances) {
		temp = temp %*% adj
		distances[which(temp > 0)] = pmin(distances[which(temp > 0)], step)
		step = step + 1	
	}

	for (i in 1:numWaterpools) {
		distances[i,i] = 0
	}

	return (distances)
}

processObservation=function(readings,probs) {
	p = rep(1, numWaterpools)
	p = p * dnorm(readings[1], mean=probs$salinity[,1], sd=(probs$salinity[,2]))
	p = p * dnorm(readings[2], mean=probs$phosphate[,1], sd=(probs$phosphate[,2]))
	p = p * dnorm(readings[3], mean=probs$nitrogen[,1], sd=(probs$nitrogen[,2]))
	return (p)
}

hmmWC=function(moveInfo,readings,positions,edges,probs) {
	if (length(moveInfo$mem) == 0) {
		adjacencyMatrix = getAdjacencyMatrix(edges)
		moveInfo$mem[[1]] = initialState(positions)
		moveInfo$mem[[2]] = transitionMatrix(adjacencyMatrix)
		moveInfo$mem[[3]] = distanceMatrix(adjacencyMatrix)
	}

	# Calculate current state vector

	# if a tourist was eaten, we can determine where the crocodile is
	eatenTourists = which(positions < 0)
	if (length(eatenTourists) > 0) {
		stateVector = rep(0, numWaterpools)
		stateVector[-eatenTourists[1]] = 1
	} else {
		stateVector = moveInfo$mem[[1]]
		stateVector = stateVector %*% moveInfo$mem[[2]]
		if (!is.na(positions[1]) && positions[1] > 0) {
			stateVector[positions[1]] = 0
		}
		if (!is.na(positions[2]) && positions[2] > 0) {
			stateVector[positions[2]] = 0
		}
 		stateVector = stateVector * processObservation(readings,probs)
		stateVector = stateVector / sum(stateVector)
	}
	# save state vector for the next iteration
	moveInfo$mem[[1]] = stateVector

	# Choosing the best action
	currPosition = positions[3]
	distances = moveInfo$mem[[3]]

	bestDistSoFar = Inf
	bestActionSoFar = NA

	# check positions leading one step ahead
	oneStepAhead = getOptions(currPosition,edges)
	for (nextPosition in oneStepAhead) {
		d = sum(distances[nextPosition,] * stateVector)
		if (d < bestDistSoFar) {
			bestActionSoFar = c(nextPosition, 0)
			bestDistSoFar = d
		}
	}

	# check positions leading two steps ahead
	for (posAfterOneStep in oneStepAhead) {
		twoStepsAhead = getOptions(posAfterOneStep,edges)
		for (nextPosition in twoStepsAhead) {
			d = sum(distances[nextPosition,] * stateVector)
			if (d < bestDistSoFar) {
				bestActionSoFar = c(posAfterOneStep, nextPosition)
				bestDistSoFar = d
			}
		}
	}	

	moveInfo$moves=bestActionSoFar

	return(moveInfo)
}
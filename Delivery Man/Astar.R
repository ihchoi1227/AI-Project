library(liqueueR)

# Check whether a package is present at (x,y)
isPackagePresent=function(x, y, packages) {
  numPackages = nrow(packages)
  for (i in 1:numPackages) {
    if (packages[i,5] == 0 # check that package status is 'not picked up'
        && isTRUE(packages[i,1] == x) && isTRUE(packages[i,2] == y)) {
      return (TRUE)
    }
  }
  return (FALSE)
}

# Determine whether current state is a goal state:
#
#   if no package is loaded, the goal is to pick a package;
#   if a package is loaded, the goal is to deliver it.
#
isGoalState=function(state, car, packages) {
  if (car$load == 0) {
    return (isPackagePresent(state$x, state$y, packages))
  }
  # Compare delivery destination coordinates to current position
  dx = packages[car$load, 3]
  dy = packages[car$load, 4]

  return (dx == state$x && dy == state$y)
}

# Determine which action leads from current position to successor position
getAction=function(currX, currY, nextX, nextY) {
  action = 5

  if (currY < nextY) {
    action = 8
  } else if (currY > nextY) {
    action = 2
  } else if (currX < nextX) {
    action = 6
  } else {
    action = 4
  }
  return (action)
}

# Return a list of legal actions.
# Walls may be blocking a move.
getLegalActions=function(state, dim) {
  # Get dimension of environment
  width = dim[1]
  height = dim[2]

  legalActions = c() # for debug purposes, '5' is omitted
  if (state$x > 1) {
    legalActions = c(legalActions, 4)
  }
  if (state$x < width) {
    legalActions = c(legalActions, 6)
  }
  if (state$y > 1) {
    legalActions = c(legalActions, 2)
  }
  if (state$y < height) {
    legalActions = c(legalActions, 8)
  }
  return (legalActions)
}

# Return successor state position resulting from the action
getSuccessor=function(state, action) {
  dx = 0
  dy = 0
  if (action == 4) {
    dx = -1
  } else if (action == 6) {
    dx = 1
  } else if (action == 2) {
    dy = -1
  } else if (action == 8) {
    dy = 1
  }
  return (list(x=state$x + dx, y=state$y + dy))
}

# Determine cost of performing 'action' being at 'currState'
getTransitionCost=function(currState, action, roads) {
  x = currState$x
  y = currState$y

  cost = 0
  if (action == 4) {
    cost = roads$vroads[x-1,y]
  } else if (action == 6) {
    cost = roads$vroads[x,y]
  } else if (action == 2) {
    cost = roads$hroads[x,y-1]
  } else if (action == 8) {
    cost = roads$hroads[x,y]
  }
  return (cost)
}

# Heuristic function
heuristic=function(state, roads, car, packages) {
  x = state$x 
  y = state$y
  numPackages = nrow(packages)

  manDists = rep(Inf, numPackages)

  if (car$load == 0) {
    offset = 0
  } else {
    offset = 2
  }

  for (i in 1:numPackages) {
    manDists[i] = abs(x - packages[i,1+offset]) + abs(y - packages[i,2+offset])
  }
  
  return (min(manDists))
}

# Priority function for a search node.
# Priority of each node is defined as (cost of reaching this node)
# + (heuristically estimated cost of path to goal).
# Note: minus (-) is used because liqueueR.PriorityQueue pops max first
priorFn=function(node) {
  return (-(node$cost + node$estimate))
}

# Returns state after performing action in current state
performAction=function(currState, action, roads, car, packages) {
  # Determine successor position
  succPosition = getSuccessor(currState, action)

  # Determine heuristic and transition cost for successor
  succEstimate = heuristic(succPosition, roads, car, packages)
  transitionCost = getTransitionCost(currState, action, roads)
        
  # Create successor state
  succState = list(
    x=succPosition$x,
    y=succPosition$y, 
    estimate=succEstimate, 
    cost=currState$cost + transitionCost + 1 # + 1 for a new turn
  )

  return (succState)
}

# Estimate cost of getting to goal if we perform given action from given state
estimateCost=function(startState, action, minCostGrid, roads, car, packages) {
  # Extract environment dimensions from roads
  width = ncol(roads$vroads)
  height = nrow(roads$hroads)

  # Create frontier as a priority queue
  frontier <- PriorityQueue$new(prioritise=priorFn)

  # Perform action 
  nextState = performAction(startState, action, roads, car, packages)
  minCostGrid[nextState$x, nextState$y] = nextState$cost

  # Add start node to frontier
  frontier$push(nextState)

  # Continue with A* algorithm until frontier is not empty
  while (frontier$size() > 0) {
    currState = frontier$pop()

    # Explore state only if there is no cheaper way to reach it  
    if (minCostGrid[currState$x, currState$y] >= currState$cost) {

      # If goal state is reached, return the cost of reaching the goal
      if (isGoalState(currState, car, packages)) {
        return (currState$cost)
      }

      # Add successors of current state to frontier
      actions = getLegalActions(currState, c(width, height))

      for (a in actions) {
        succState = performAction(currState, a, roads, car, packages)

        # Add state to frontier only if it was reached with a better cost than before
        if (minCostGrid[succState$x, succState$y] > succState$cost) {
          frontier$push(succState)
          minCostGrid[succState$x, succState$y] = succState$cost
        }
      }
    }
  }

  return (Inf)
}

# A* search implementation
astarDM=function(roads, car, packages) {
  # Extract environment dimensions from roads
  width = ncol(roads$vroads)
  height = nrow(roads$hroads)

  # Grid to save minimal cost of getting to a node
  # Used to avoid expanding extra states
  minCostGrid = matrix(Inf, nrow=height, ncol=width)

  # Create a node for start state
  startState = list(
    x=car$x, 
    y=car$y,
    estimate=0, 
    cost=0 
  )
  startEstimate = heuristic(startState, roads, car, packages)
  startState$estimate = startEstimate
  
  # Mark cost of getting to start state as 0
  minCostGrid[startState$x,startState$y] = 0

  actions = getLegalActions(startState, c(width, height))

  bestActionSoFar = 5
  bestCostSoFar = Inf # stop is never a rational choice

  for (a in actions) {
    currCost = estimateCost(startState, a, minCostGrid, roads, car, packages)
    if (currCost < bestCostSoFar) {
      bestCostSoFar = currCost
      bestActionSoFar = a
    }
  }

  car$nextMove = bestActionSoFar
  return (car)
}
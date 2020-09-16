##~~Patch Notes
#Fix drop units errors +
#Forager, brwForager Paths as Linestrings +
#Allow multiple linestrings in Forager paths +
#Add crwForager -
#dont relocate forager in ArrayEnvironment after final trial +

#' @import sp
#' @import dplyr
#' @importFrom circular circular rwrappedcauchy
#' @importFrom purrr reduce
#' @importFrom stats runif
#' @importFrom units set_units
#' @import methods
#' @import sf



#~~~~~~~~~~~~~~~~~~~~~~Object Set-up~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Environment Constructor
#' @description  Constructor Function for objects of Reference Class 'Environment'
#' @field foragers a list containing Reference Class objects that inherit from 'Forager' class.
#' @field patches a simple feature data frame with geometries representing patches in the geometry column
#' @field bounds a spatial features collection with one or more polygons defining the region in which any foragers can be created and moved to
#' @export
#' @export Environment
#' @examples
#' \dontrun{
#' Environment()
#' newForagers <- list(Forager(), Forager())
#' #foragers are created by default between -50 and 50 on both the x and y axis
#' newPatches <- lapply(rep(2, times = 20), function(x) st_point(runif(x, min = -50, max = 50)))
#' newPatches <- newPatches %>% st_sfc %>% data.frame(geom = ., NAME = as.character(1:length(.)))
#' newPatches <- newPatches %>% st_sf() %>% st_buffer(1)
#' newBounds <- st_sfc(st_convex_hull(reduce(newPatches$geometry, c)), crs = st_crs(newPatches))
#' Environment(newForagers, newPatches, newBounds)
#' }
#' @import methods
#' @import dplyr
Environment <- setRefClass("Environment",
                           field = list(foragers = "list", patches = "sf", bounds = "sfc_POLYGON"),
                           method = list(
                             initialize = function(...,
                                                   patches = lapply(rep(2, times = 10), function(x) st_point(runif(x, min = -15, max = 15))) %>% st_sfc %>% data.frame(geom = ., NAME = as.character(1:length(.))) %>% st_sf() %>% st_buffer(1) %>% st_set_crs(32610),
                                                   bounds = st_sfc(st_buffer(st_convex_hull(reduce(patches$geometry, c)),1), crs = st_crs(patches)), #default bounds created as convex hull around patches with a buffer of one
                                                   foragers = createForagers(3, bounds = bounds, speed = 2, quiet = TRUE) #default creation of forager within bounds
                             ){  "Set default values for variables that are not entered manually"
                               callSuper(..., foragers = foragers, patches = patches, bounds = bounds)
                             },
                             progress = function(){ "moves the entire environment forward one timestep by moving all foragers once and processing impact of foraging on patches"
                               for (forager in foragers[sample(length(foragers))]) {#operates on each forager in random order
                                 if (! forager$targeting) forager$setTarget(patches) #If the forager does not currently have a target, look for one
                                 forager$move(bounds = bounds)
                               }
                             },
                             plotPatches = function(){ "displays the current patch geometries"
                               return(ggplot(data = patches) + geom_sf(fill = "green", axes = TRUE) + theme_classic()) #plots all geometries in patches field
                             },
                             plotCurrent = function(){ "displays the current location of all foragers and patch geometries"
                               fLocations <- st_sf(geom = reduce(lapply(foragers, function(X) X$location),c), TYPE = as.character(sapply(foragers, class)), SPEED = sapply(foragers, function(X) X$speed)) #gets simple feature geometry column for each forager location in a list, then combines all into one simple feature geometry column
                               cPlot <- ggplot() + geom_sf(data = patches, fill = "green") + geom_sf(data = fLocations, aes(color = TYPE), shape = 13, size = 3) + geom_sf(data = bounds, linetype = "longdash", fill = NA) + theme_classic()
                               return(cPlot)
                             },
                             plotPaths = function(){ "displays the current location and path of all foragers as well as patch geometries"
                               cPlot <- plotCurrent()
                               paths <- purrr::reduce(lapply(foragers, function(X) st_cast(X$path, "LINESTRING")), c)
                               pathsPlot <- cPlot + geom_sf(data = paths, color = 1:length(paths))
                               return(pathsPlot)
                             })
)

#' @title ArrayEnvironment Constructor
#
#' @description Constructor function for the ArrayEnvironment class, which inherits from the Environment class. In addition to the features of a normal environment, this object will reset foragers once all patches have been visited (a trial), while tracking a single foragers sequences over all trials
#' @inherit Environment
#' @field foragers a list containing Reference Class objects that inherit from 'Forager' class.
#' @field patches a simple feature data frame with geometries representing patches in the geometry column
#' @field bounds a spatial features collection with one or more polygons defining the region in which any foragers can be created and moved to
#' @field sequences the names of all patches visited (in order) by foragers (only really works with a single forager currently, sequences by multiple foragers will be mixed together). Usually the defualt of an empty character is used
#' @field array the name of the array being run, if it has one
#' @field trials the number of trials that have been run. Has no effect on the environment, merely usefull as metadata. Typically uses default value of 0 at creation.
#' @export
#' @export ArrayEnvironment
#' @import methods
ArrayEnvironment <- setRefClass("ArrayEnvironment", fields= list(sequence = "character", array = "character", trials = "numeric"), contains= "Environment",
                                method = list(
                                  initialize = function(...,
                                                        sequence = character(0),
                                                        array = "DT",
                                                        trials = 0) { "Set default values for variables that are not entered manually"
                                    callSuper(..., sequence = sequence, array = array, trials = trials)
                                    if (length(foragers) > 1) {
                                      warning("ArrayEnvironment can only handle a single forager. All additional foragers removed")
                                      foragers <<- foragers[1]
                                    }
                                  },
                                  progress = function(){ "Moves all foragers as desribed in parent class Environment. Then looks to see if any foragers have either visited all patches or moved 10000 steps in their current trial.
                                    If any forager has, it's patch visitation sequence is combined with the sequence field, its location field is randomized within bounds, and its visitSeq field is reset"
                                    if (length(foragers[[1]]$path) < trials) {
                                      foragers[[1]]$location <<- generateBoundedPoint(bounds)
                                      foragers[[1]]$path[[trials]] <<- st_cast(test$foragers[[1]]$location, "LINESTRING")[[1]]
                                      foragers[[1]]$bearing <<- as.numeric(rwrappedcauchy(n = 1, mu = circular(0), rho = 0))
                                      foragers[[1]]$visitSeq <<- rep("NA", times = foragers[[1]]$repeatAvoid)
                                    }
                                    callSuper()
                                    for (forager in foragers) {
                                      if (length(unique(forager$visitSeq[-c(1:forager$repeatAvoid)])) == nrow(patches) | length(forager$path[[1]]) >= 5000) { #end conditions for trial 1) Forager has been to all patches, 2) forager has made 2500 steps (length of multipoints object in path[[1]] counts both x and y coordinates, so length of 5000 is == 2500 steps)
                                        sequence <<- c(sequence, forager$visitSeq[-c(1:forager$repeatAvoid)])
                                        trials <<- trials + 1
                                      }
                                    }
                                  }
                                ))

#' @title Forager Constructor
#
#' @description Constructor function for Forager class objects. Objects of this class have methods for finding targets, moving, and plotting. Movement switches between brownian motion (with step lengths drawn from a gamme distribution) and directed toward a target
#' @field location a simple features collection with a single POINT class object
#' @field bearing a numeric that gives the current bearing in radians. Default starting bearing is drawn from random uniform circular distribution
#' @field speed a numeric value that gives the scale parameter of the gamma distribution from which step lengths are draw. Because the shape parameter of this distribution is set to 1, the speed variable will equal the average step length
#' @field sight a numeric giving the distance at which the forager object can detect patches when in an environment
#' @field path a simple features collection with a single multipoint object containing previous locations of the forager. Usually objects are created with default multipoint objects containing only the initial location
#' @field visitSeq a character vector with the names (in order) of all patches visited. Must start with NAs equal in number to to the repeatAvoid variable (this is done by default if no value is given)
#' @field targeting a logical giving which mode of movement the forager is in. Initial value usually uses default "FALSE"
#' @field repeatAvoid numeric giving the number of different patches a forager must visit before targeting a recently visited patch again
#' @field target a simple feature data frame with a single row containing the target patch location and values
#' @export
#' @export Forager
#' @import methods
#' @importFrom circular circular rcircularuniform
#' @import sp
#' @import ggplot2
Forager <- setRefClass("Forager",
                       field = list(location="sfc", bearing ="numeric", speed = "numeric", sight = "numeric", path = "sfc", visitSeq = "character", targeting = "logical", repeatAvoid = "numeric", target = "sf"),
                       method = list(initialize = function(...,location = st_sfc(st_point(runif(2, -50, 50))),
                                                           bearing = as.numeric(rcircularuniform(1)),
                                                           speed = 1, #the mean distance traveled when moving randomly, drawn from a gamma distribution with k (shape) = 1, and theta(scale) = speed. When targeting, moves at 2X speed
                                                           sight = 4,
                                                           path = st_cast(location, "LINESTRING"),
                                                           repeatAvoid = 2,
                                                           visitSeq = rep("NA", times = repeatAvoid),#NA's prevent errors when checking for recent visits
                                                           targeting = FALSE,
                                                           target = st_sf(st_sfc(st_multipolygon()))) {
                         callSuper(..., location = location, bearing = bearing, speed = speed, sight = sight, path = path, visitSeq = visitSeq, targeting = targeting, target = target, repeatAvoid = repeatAvoid)
                       },
                       setTarget = function(patches) { "checks if the location of any simple features in the geometry column of the patches argument (class sfc) are within sight range.
                         If so, sets target to a patch randomly selected from those within sight and then sets the targetting field to TRUE"
                         if (! "sf" %in% class(patches)) stop("the 'patches' argument of setTarget must be a simple feature collection")
                         choices <- getChoices(.self, patches)
                         if (! "sf" %in% class(choices)) return()
                         target <<- tSelect(choices)
                         targeting <<- TRUE
                       },
                       tSelect = function(choices) { "calculates probablilites and selects a target from a list of options but DOES NOT set the target (use setTarget for this, which calls tSelect). A forager will always select a patch it is on if possible. Otherwise, selection is random"
                         if (sum(choices$DIST == 0) > 0) return(choices[which(choices$DIST == 0),][1,])     #if on a patch, return that patch as target. If on multiple, return the first
                         choices$prob <- rep(1/nrow(choices), times = nrow(choices)) #random choice
                         return(selector(choices))
                       },
                       move = function(bounds = NA) {"moves the forager using rules given in class description. Updates relevant fields including location, bearing, visitSeq, bearing, targetting, and path"
                         if (targeting) { #if a target is given, assess proximity and action
                           bearing <<- lineBearing(st_nearest_points(st_set_crs(st_sfc(location), st_crs(target)), target)) #set bearing toward target
                           if (abs(set_units(st_distance(location, target$geometry), NULL)) <= 2 * speed) { #if target is in reach
                             location[1] <<- st_cast(st_nearest_points(location, target$geometry), "POINT")[[2]] # set location to patch location
                             visitSeq <<- c(visitSeq, as.character(target$NAME)) #add patch to visitSeq
                             targeting <<- FALSE
                           } else location[1] <<- location + 2 * speed * c(cos(bearing), sin(bearing)) #otherwise, move closer
                         } else { #if not targeting
                           bearing <<- as.numeric(rcircularuniform(1))
                           location[1] <<- location + rgamma(1, shape = 1, scale = speed) * c(cos(bearing), sin(bearing))
                           if (!is.na(bounds)){ #check for valid bounds
                             while(! st_within(location, bounds, sparse = FALSE)[1,1]) {#check if out of bounds
                               location[1] <<- st_point(path[[length(path)]][nrow(path[[length(path)]]),]) #reset location
                               bearing <<- as.numeric(rcircularuniform(1))
                               location[1] <<- location + rgamma(1, shape = 1, scale = speed) * c(cos(bearing), sin(bearing))
                             }
                           }
                         }
                         path[[length(path)]] <<- st_linestring(rbind(path[[length(path)]], location[[1]]))
                         attr(path, "bbox") <<- setNames(unlist(bbox(as_Spatial(path))), c("xmin", "ymin", "xmax", "ymax"))
                         attr(location, "bbox") <<- setNames(unlist(bbox(as_Spatial(location))), c("xmin", "ymin", "xmax", "ymax"))
                       },
                       plot = function(){ "plots the location and path of the forager"
                         fPlot <-ggplot() + geom_sf(data = location, shape = 13, size = 4) + geom_sf(data = st_cast(path, "LINESTRING")) + theme_classic()
                         return(fPlot)
                       })
)

#' @title brwForager constructor
#
#' @description Constructor function for a Biased Random Walk Forager class objects. Objects of this class have methods for finding targets, moving, and plotting. Movement switches between a biased random walk (with step lengths drawn from a gamma distribution and bearing drawn from a wrapped cauchy) and directed motion toward a target
#' @inherit Forager
#' @field location a simple features collection with a single POINT class object
#' @field bearing a numeric that gives the current bearing in radians. Default starting bearing is drawn from random uniform circular distribution
#' @field speed a numeric value that gives the scale parameter of the gamma distribution from which step lengths are draw. Because the shape parameter of this distribution is set to 1, the speed variable will equal the average step length
#' @field sight a numeric giving the distance at which the forager object can detect patches when in an environment
#' @field path a simple features collection with a single multipoint object containing previous locations of the forager. Usually objects are created with default multipoint objects containing only the initial location
#' @field visitSeq a character vector with the names (in order) of all patches visited. Must start with NAs equal in number to to the repeatAvoid variable (this is done by default if no value is given)
#' @field targeting a logical giving which mode of movement the forager is in. Initial value usually uses default "FALSE"
#' @field repeatAvoid numeric giving the number of different patches a forager must visit before targeting a recently visited patch again
#' @field target a simple feature data frame with a single row containing the target patch location and values
#' @field turnBias a numeric giving the radians from which the center of the wrapped cauchy distribution from which new bearings are drawn should be shifted from the previous bearing.
#' @field persistence a numeric between 0 and one which gives the concentration of the wrapped cauchy distribution from which new bearings are drawn, where a 1 results in a turning angle equal to turnBias every step and a 0 results in a circular random uniform probability distribution
#' @export
#' @export BRWForager
#' @import methods
#' @importFrom circular circular rcircularuniform rwrappedcauchy
#' @import sp
#' @import ggplot2
BRWForager <- setRefClass("brwForager", fields = list(turnBias = "numeric", persistence = "numeric"), contains = "Forager",
                          methods = list(
                            initialize = function(..., turnBias = 0, persistence = 0.6) {
                              callSuper(..., turnBias = turnBias, persistence = persistence)
                            },
                            move = function(bounds = NA){
                              if (targeting) { #if a target is given, assess proximity and action
                                bearing <<- lineBearing(st_nearest_points(st_set_crs(st_sfc(location), st_crs(target)), target)) #set bearing toward target
                                if (abs(set_units(st_distance(location, target$geometry), NULL)) <= 2 * speed) { #if target is in reach
                                  location[1] <<- st_cast(st_nearest_points(location, target$geometry), "POINT")[[2]] # set location to patch location
                                  visitSeq <<- c(visitSeq, as.character(target$NAME)) #add patch to visitSeq
                                  targeting <<- FALSE
                                } else location[1] <<- location + speed * c(cos(bearing), sin(bearing)) #otherwise, move closer
                              } else { #if not targeting
                                bearing <<- as.numeric(circular(bearing) + rwrappedcauchy(n = 1, mu = as.circular(turnBias, type = "angles", units = "radians", template = "none", zero = 0, modulo = "asis", rotation = "counter"), rho = persistence))
                                location[1] <<- location + speed * c(cos(bearing), sin(bearing))
                                if (!is.na(bounds)){ #check for valid bounds
                                  turnVarIncrease <- 0
                                  while(! st_within(location, bounds, sparse = FALSE)[1,1]) {#check if out of bounds
                                    location[1] <<- st_point(path[[length(path)]][nrow(path[[length(path)]]),]) #reset location
                                    if(persistence-turnVarIncrease > 0.2) turnVarIncrease <- turnVarIncrease + 0.02 #relax directional persistence
                                    bearing <<- as.numeric(circular(bearing) + rwrappedcauchy(n = 1, mu = circular(0), rho = persistence - turnVarIncrease)) # get new bearing
                                    location[1] <<- location + speed * c(cos(bearing), sin(bearing)) #try moving again
                                  }
                                }
                              }
                              path[[length(path)]] <<- st_linestring(rbind(path[[length(path)]], location[[1]]))
                              attr(path, "bbox") <<- setNames(unlist(bbox(as_Spatial(path))), c("xmin", "ymin", "xmax", "ymax"))
                              attr(location, "bbox") <<- setNames(unlist(bbox(as_Spatial(location))), c("xmin", "ymin", "xmax", "ymax"))
                            })
)

#' @title dForager Constrctor
#
#' @description Constructor function for distance Forager class objects. Obects of this class behave exactly like forager class objects, but when selecting a target, distance foragers prefer closer targets, discounting targets probability of selection by the cube of their distance
#' @inherit Forager
#' @field location A simple features collection with a single POINT class object
#' @field bearing A numeric that gives the current bearing in radians. Default starting bearing is drawn from random uniform circular distribution
#' @field speed A numeric value that gives the scale parameter of the gamma distribution from which step lengths are draw. Because the shape parameter of this distribution is set to 1, the speed variable will equal the average step length
#' @field sight A numeric giving the distance at which the forager object can detect patches when in an environment
#' @field path A simple features collection with a single multipoint object containing previous locations of the forager. Usually objects are created with default multipoint objects containing only the initial location
#' @field visitSeq A character vector with the names (in order) of all patches visited. Must start with NAs equal in number to to the repeatAvoid variable (this is done by default if no value is given)
#' @field targeting A logical giving which mode of movement the forager is in. Initial value usually uses default "FALSE"
#' @field repeatAvoid Numeric giving the number of different patches a forager must visit before targeting a recently visited patch again
#' @field target A simple feature data frame with a single row containing the target patch location and values
#' @export
#' @export dForager
#' @import methods
dForager <- setRefClass("distanceForager", fields = list(), contains = "Forager",
                        methods = list(
                          tSelect = function(choices) {
                            if (sum(choices$DIST == 0) > 0) return(choices[which(choices$DIST == 0),][1,])                                         #if on a patch, return that patch as target. If on multiple, return the first
                            choices$prob <- choices$DIST^3/sum(choices$DIST^3) #distance discounted choice
                            return(selector(choices))
                          }))

#' @title ddbrwForager constructor
#
#' @description Constructor function for a Distance Discounting Biased Random Walk Forager class objects. Objects of this class have methods for finding targets, moving, and plotting. Movement switches between a biased random walk (with step lengths drawn from a gamma distribution and bearing drawn from a wrapped cauchy) and directed motion toward a target. Target selection is wieghted toward closer targets.
#' @inherit brwForager
#' @field location a simple features collection with a single POINT class object
#' @field bearing a numeric that gives the current bearing in radians. Default starting bearing is drawn from random uniform circular distribution
#' @field speed a numeric value that gives the scale parameter of the gamma distribution from which step lengths are draw. Because the shape parameter of this distribution is set to 1, the speed variable will equal the average step length
#' @field sight a numeric giving the distance at which the forager object can detect patches when in an environment
#' @field path a simple features collection with a single multipoint object containing previous locations of the forager. Usually objects are created with default multipoint objects containing only the initial location
#' @field visitSeq a character vector with the names (in order) of all patches visited. Must start with NAs equal in number to to the repeatAvoid variable (this is done by default if no value is given)
#' @field targeting a logical giving which mode of movement the forager is in. Initial value usually uses default "FALSE"
#' @field repeatAvoid numeric giving the number of different patches a forager must visit before targeting a recently visited patch again
#' @field target a simple feature data frame with a single row containing the target patch location and values
#' @field turnBias a numeric giving the radians from which the center of the wrapped cauchy distribution from which new bearings are drawn should be shifted from the previous bearing.
#' @field persistence a numeric between 0 and one which gives the concentration of the wrapped cauchy distribution from which new bearings are drawn, where a 1 results in a turning angle equal to turnBias every step and a 0 results in a circular random uniform probability distribution
#' @export
#' @export DDBRWForager
#' @import methods
#' @importFrom circular circular rcircularuniform rwrappedcauchy
#' @import sp
#' @import ggplot2
DDBRWForager <- setRefClass("ddbrwForager", contains = "brwForager",
                          methods = list(
                            initialize = function(...) {
                              callSuper(...)
                            },
                            tSelect = function(choices) {
                              if (sum(choices$DIST == 0) > 0) return(choices[which(choices$DIST == 0),][1,])                                         #if on a patch, return that patch as target. If on multiple, return the first
                              choices$prob <- choices$DIST^3/sum(choices$DIST^3) #distance discounted choice
                              return(selector(choices))
                            })
)
# Function to make a selection given a set of probabilities
#
#' @param choices a datafame with a column named "prob". Values in this column should be numeric and and up to 1.
selector <- function(choices){
  randomizer <- runif(1)
  i = 1
  while(choices$prob[i] <= randomizer) {
    randomizer <- randomizer - choices$prob[i]
    i <- i + 1
  }
  return(choices[i,])
}

# Function to calculate distance between a forager and a set of patches, then return the patches withing sight range of the forager and their distances
#
#' @param forager A single object of Forager class
#' @param patches A simple features data frame with a geom column containing the geometries of patches.
#' @importFrom units set_units
getChoices <- function(forager, patches) {
  recentVisits <- forager$visitSeq[c((length(forager$visitSeq) - (forager$repeatAvoid - 1)):length(forager$visitSeq))] #check which patches forager has visited recently
  patches <- patches[! patches$NAME %in% recentVisits,] #remove recently visited patches
  distances <- set_units(st_distance(forager$location, patches), NULL) #get distances to each patch
  if (sum(distances <= forager$sight) == 0) return(NA)                                                                        #if no patches in sight, return no target
  choices <- patches[which(distances[1,] <= forager$sight),]
  choices$DIST <- distances[distances[1,] <= forager$sight]
  return(choices)
}

# Function to get the bearing from the start point of a linestring to the endpoint
#
#' @param linestring A st_LINESTRING object
lineBearing <- function(linestring) {
  endpoints <- st_cast(linestring, "POINT")
  deltaXY <- endpoints[[2]] - endpoints[[1]]
  return(atan2(x = deltaXY[1], y = deltaXY[2]))
}

#' @title Forager Constructor Wrapper
#' @description A wrapper function to create a list of foragers from lists of parameter values. Handles single values if parameter is the same for all foragers. Assigns default values for empty fields.
#
#' @param numForagers A numeric giving the number of foragers to be created
#' @param type A character giving the type of forager to create. Currenly handles "Random" (default, parent class Forager) and "BRW" (brwForager)
#' @param bounds A simple features collection of polygons within which all foragers must be created
#' @param locations A list of simple features collections each containing a single st_POINT giving the initial location of a new forager. If NA or empty, foragers will be created at random locations within bounds
#' @param bearings A single numeric giving the initial bearing of all new foragers in radians or a list of numerics equal in length to the value of numForagers
#' @param speeds A single numeric giving the initial speed of all new foragers or a list of numerics equal in length to the value of numForagers
#' @param persistences A single numeric giving the initial persistence of all new foragers or a list of numerics equal in length to the value of numForagers. Ignored unless type = "BRW"
#' @param sights A single numeric giving the initial sight range of all new foragers or a list of numerics equal in length to the value of numForagers
#' @param repeatAvoids A single numeric giving the number of unique patches a forager must visit before revisting a patch, or a list of numerics equal in length to the value of numForagers
#' @param quiet A logical indicating whether to suppress warnings about default values used to fill in empty fields
#' @param CRS A numeric giving the CRS code to assign to forager locations. If bounds argument is defined, the CRS of that object will be used and this argument will be ignored
#' @param turnBiases A single numeric giving the average turn angle of a forager in radians, or a list of numerics equal in length to the value of numForagers. Ignored unless type = "BRW"
#' @export
#' @importFrom circular circular rcircularuniform rwrappedcauchy
createForagers <- function(numForagers, type = "Random", bounds = NA, locations = NA, bearings = NA, speeds = NA, persistences = NA, sights = NA, repeatAvoids = NA, quiet = FALSE, CRS = "NA", turnBiases = NA) {
  if (is.na(bearings)) {
    if (! quiet) warning("No bearings given, initial values drawn from circular random uniform distribution")
    bearings <- as.numeric(rwrappedcauchy(n = numForagers, mu = circular(0), rho = 0))
  }
  if(is.na(speeds)) {
    if (! quiet) warning("No speeds given, initial values set to 1")
    speeds <- 1
  }
  if (is.na(sights)) {
    if (! quiet) warning("No sight ranges given, initial values set to 5")
    sights <- 5
  }
  if (is.na(repeatAvoids)) {
    if (! quiet) warning("No patch avoidance memory given. Initial values set to avoid 2 most recently visited patches")
    repeatAvoids <-  2
  }
  if (is.na(locations)) {
    if (is.na(bounds)) {
      if (! quiet) warning("No bounds or locations given. Bounds set from -50 to 50 for x and y axes by default")
      bounds <- st_sfc(st_polygon(list(matrix(c(50, 50, -50, 50, -50, -50, 50, -50, 50, 50), ncol = 2, byrow = TRUE))), crs = CRS)
    }
    locations <- vector("list", length = numForagers)
    for(i in 1:length(locations)) locations[[i]] <- generateBoundedPoint(bounds)
  } else {
    if ("data.frame" %in% class(locations)) {
      if(!ncol(locations) == 2) stop("locations argument must be either a list of coordinate pairs, a list of spatial points, or a dataframe with 2 columns for x and y coodinates respectively")
      locationPoints <- vector("list", length = nrow(locations))
      for (i in 1:length(locationPoints)) locationsPoints[[i]] <- st_point(locations[i,])
      locations <- locationPoints
    }
  }
  if (type == "BRW"){
     if(is.na(persistences)){
       if (! quiet) warning("No turning persistences given, initial values set to 0.7")
       persistences <- 0.7
     }
    if(is.na(turnBiases)){
      if (! quiet) warning("No turning biases given, initial values set to 0")
      turnBiases <- 0
    }
  }
  if(!length(locations) == numForagers) stop("number of locations given must equal value of numForagers argument")
  locations <- lapply(locations, to_sf_point, crs = st_crs(bounds))
  foragers <- vector("list", length = numForagers)
  parameters <- st_sf(geom = reduce(locations,c), BEARING = bearings, SPEED = speeds, PERSISTENCE = persistences, SIGHT = sights, REPEATAVOID = repeatAvoids, BIAS = turnBiases)
  if(type == "Random") for (i in 1:numForagers) foragers[[i]] <- Forager(location=parameters$geom[i], bearing = parameters$BEARING[i], speed = parameters$SPEED[i], sight = parameters$SIGHT[i], repeatAvoid = parameters$REPEATAVOID[i])
  if(type == "BRW") for (i in 1:numForagers) foragers[[i]] <- BRWForager(location=parameters$geom[i], bearing = parameters$BEARING[i], speed = parameters$SPEED[i], sight = parameters$SIGHT[i], repeatAvoid = parameters$REPEATAVOID[i], persistence = parameters$PERSISTENCE[i], turnBias = parameters$BIAS[i])
  return(foragers)
}

# Function to flexibly handle multiple point formats and converts to a simple features collection
#
#' @param point A length two numeric vector, an st_POINT, or and sfc_POINT class object
#' @param crs A numeric giving the CRS code to apply to the point argument
to_sf_point <- function(point, crs = NA) {
  if (is.numeric(point)) {
    if (!length(point) == 2) stop("point argument is numeric, but not of length 2")
    point <- st_point(point)
  }
  if ("POINT" %in% class(point)) point <- st_sfc(point, crs = crs)
  if ("sfc_POINT" %in% class(point)) return(point) else stop("point argument must be a numeric vector, a POINT object, or a geometry collection with POINT objects")
}


# Function to create a simple features collection containing a single point created randomly within a set of polygons
#
#' @param bounds A simple features collection containing polygons within which new point must be created.
#' @importFrom stats runif
generateBoundedPoint <- function(bounds) {
  location <- c(runif(1, st_bbox(bounds)["xmin"], st_bbox(bounds)["xmax"]), runif(1, st_bbox(bounds)["ymin"], st_bbox(bounds)["ymax"])) %>% st_point() %>% st_sfc(crs = st_crs(bounds))
  while (!st_within(location, bounds, sparse = FALSE)[1,1]) location <- c(runif(1, st_bbox(bounds)["xmin"], st_bbox(bounds)["xmax"]), runif(1, st_bbox(bounds)["ymin"], st_bbox(bounds)["ymax"])) %>% st_point() %>% st_sfc(crs = st_crs(bounds)) #if point not within bounds, recreate
  return(location)
}


###~~~~~~Sample Script~~~~~~~~

#' @examples
#' newPatches <- lapply(rep(2, times = 20), function(x) st_point(runif(x, min = -100, max = 100))) %>% st_sfc %>% data.frame(geom = ., NAME = as.character(1:length(.))) %>% st_sf() %>% st_buffer(1) %>% st_set_crs(32610) #creates sf dataframe with 200 randomly created circular patches with radius of 2
#' brwForagers <- createForagers(5, type = "BRW", bounds = st_sfc(st_convex_hull(reduce(testPatches$geometry, c)), crs = st_crs(testPatches)), speeds = 4, sights = 8, quiet = TRUE)
#' rForagers <- createForagers(5, bounds = st_sfc(st_convex_hull(reduce(testPatches$geometry, c)), crs = st_crs(testPatches)), speeds = c(2,4,6,8,10), sights = 5, quiet = TRUE)
#' mixedForagers <- c(brwForagers, rForagers)
#' testEnviron <- Environment(foragers = testForagers, patches = testPatches)
#' testEnviron$plotCurrent()
#' for (t in seq(50)) testEnviron$progress()
#' testEnviron$plotPaths()


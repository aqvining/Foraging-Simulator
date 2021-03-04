#' @title Forager Constructor
#
#' @description Constructor function for Forager class objects. Objects of this class have methods for finding targets, moving, and plotting. Movement switches between brownian motion (with step lengths drawn from a gamme distribution) and directed toward a target
#' @field name a character
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
                       field = list(name = "character",
                                    location="sfc",
                                    bearing ="numeric",
                                    speed = "numeric",
                                    sight = "numeric",
                                    path = "sfc",
                                    visitSeq = "character",
                                    targeting = "logical",
                                    foraging = "logical",
                                    giving_up_density = "numeric",
                                    efficiency = "numeric",
                                    energy = "numeric",
                                    repeatAvoid = "numeric",
                                    choice_determinism = "numeric",
                                    target = "sf"),
                       method = list(
                         initialize = function(...,
                                               name = "Default Forager",
                                               location = st_sfc(st_point(runif(2, -50, 50))),
                                               bearing = as.numeric(rcircularuniform(1)),
                                               speed = 1, #the mean distance traveled when moving randomly, drawn from a gamma distribution with k (shape) = 1, and theta(scale) = speed. When targeting, moves at 2X speed
                                               sight = 4,
                                               path = st_cast(location, "LINESTRING"),
                                               visitSeq = rep("NA", times = repeatAvoid),#NA's prevent errors when checking for recent visits
                                               targeting = FALSE,
                                               foraging = FALSE,
                                               giving_up_density = 0.5,
                                               efficiency = 0.1,
                                               energy = 0,
                                               repeatAvoid = 0,
                                               choice_determinism = 0,
                                               target = st_sf(st_sfc(st_multipolygon()))) {
                           callSuper(...,
                                     name = name,
                                     location = location,
                                     bearing = bearing,
                                     speed = speed,
                                     sight = sight,
                                     path = path,
                                     visitSeq = visitSeq,
                                     targeting = targeting,
                                     foraging = foraging,
                                     giving_up_density = giving_up_density,
                                     efficiency = efficiency,
                                     energy = energy,
                                     target = target,
                                     repeatAvoid = repeatAvoid,
                                     choice_determinism = choice_determinism)
                           if (speed < 0) stop("speed must be a postive value")
                         },

                         setTarget = function(patches) { "checks if the location of any simple features in the geometry column of the patches argument (class sfc) are within sight range.
                         If so, sets target to a patch randomly selected from those within sight and then sets the targetting field to TRUE"
                           if (! "sf" %in% class(patches)) stop("the 'patches' argument of setTarget must be a simple feature collection")
                           choices <- getChoices(.self, patches)
                           if (! "sf" %in% class(choices)) return() #no target set if no patches found in sight
                           target <<- tSelect(choices)
                           targeting <<- TRUE
                         },

                         tSelect = function(choices) { "calculates probablilites and selects a target from a list of options but DOES NOT set the target (use setTarget for this, which calls tSelect)."
                           #Choices is a spatial features df with same structure as patches plus a DIST column giving the patch distance to forager.
                           choices$attraction <- (choices$VALUE)/(choices$DIST + 0.001) #Sets attraction based on value and distance. Offset is to avoid dividing by 0 if on a patch boundary
                           choices$attraction <- scale_attraction(choices$attraction, choice_determinism) #scale attractions based on choice determinism
                           choices$prob <- choices$attraction / sum(choices$attraction)
                           return(selector(choices))
                         },

                         move = function(bounds = NA) {"moves the forager using rules given in class description. Updates relevant fields including location, bearing, visitSeq, bearing, targetting, and path"
                           if (targeting) { #if a target is given, assess proximity and action
                             bearing <<- lineBearing(st_nearest_points(st_set_crs(st_sfc(location), st_crs(target)), target)) #set bearing toward target
                             if (abs(set_units(st_distance(location, target$geometry), NULL)) <= 2 * speed) { #if target is in reach
                               location[1] <<- st_cast(st_nearest_points(location, target$geometry), "POINT")[[2]] # set location to patch location
                               visitSeq <<- c(visitSeq, as.character(target$NAME)) #add patch to visitSeq
                               targeting <<- FALSE
                               foraging <<- TRUE
                             } else location[1] <<- location + 2 * speed * c(cos(bearing), sin(bearing)) #otherwise, move closer

                           } else { #if not targeting
                             if(foraging) bounds = target$geometry
                             bearing <<- as.numeric(rcircularuniform(1))
                             location[1] <<- location + rgamma(1,
                                                               shape = 1,
                                                               scale = speed/(1 + foraging)) * c(cos(bearing), sin(bearing))
                             if (!is.na(bounds)){ #check for valid bounds
                               while(! st_within(location, bounds, sparse = FALSE)[1,1]) {#check if out of bounds
                                 location[1] <<- st_point(path[[length(path)]][nrow(path[[length(path)]]),]) #reset location: gets the last linestring in the path sfc, then gets the last row of that linestring matrix
                                 bearing <<- as.numeric(rcircularuniform(1))
                                 location[1] <<- location + rgamma(1,
                                                                   shape = 1,
                                                                   scale = speed/(1 + foraging)) * c(cos(bearing), sin(bearing)) #speed reduced by half if foraging
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
#' @description Constructor function for a Biased Random Walk Forager class objects. Objects of this class have methods for finding targets, moving, and plotting. Movement switches between a biased random walk (with step lengths drawn from a gamma distribution and bearing drawn from a wrapped cauchy) and directed motion toward a target (retaining concentration, but not bias)
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
#' @field concentration a numeric between 0 and one which gives the concentration of the wrapped cauchy distribution from which new bearings are drawn, where a 1 results in a turning angle equal to turnBias every step and a 0 results in a circular random uniform probability distribution
#' @field persistence a numeric between 0 and 1 which weights the value of last step direction relative to the direction of a target when selecting new bearing
#' @export
#' @export BRWForager
#' @import methods
#' @importFrom circular circular rcircularuniform rwrappedcauchy
#' @import sp
#' @import ggplot2
BRWForager <- setRefClass("brwForager", fields = list(turnBias = "numeric", concentration = "numeric", persistence = "numeric"), contains = "Forager",
                          methods = list(
                            initialize = function(..., turnBias = 0, concentration = 0.6, persistence = 0.5) {
                              callSuper(..., turnBias = turnBias, concentration = concentration, persistence = persistence)
                            },
                            move = function(bounds = NA){
                              #This function gets new x, y position for the forager.
                              if (foraging) {
                                callSuper()
                                return()
                              }

                              if (targeting) { #if a target is given, assess proximity and action
                                target_direction <- location %>% st_set_crs(value = st_crs(target)) %>% #prepares location for comparison to target (which may have multiple location points)
                                  st_nearest_points(target) %>% #returns linestring between forager and target
                                  lineBearing() #gets angle between endpoints of linestring

                                if (abs(set_units(st_distance(location, target$geometry), NULL)) <= 2 * speed) { #if target is in reach
                                  location[1] <<- st_cast(st_nearest_points(location, target$geometry), "POINT")[[2]] # set location to patch location
                                  visitSeq <<- c(visitSeq, as.character(target$NAME)) #add patch to visitSeq
                                  bearing <<- target_direction
                                  targeting <<- FALSE
                                  foraging <<- TRUE

                                } else { #if target is not in reach
                                  target_deviation <- circular(target_direction - bearing, modulo = "2pi")
                                  if(abs(target_deviation) > pi) target_deviation <- target_deviation - (2*pi)
                                  targeting_bias <- (1-persistence) * target_deviation #peak of turning angle probability density should be between current bearing and target direction, with location in this range determined by persistence
                                  bearing <<- as.numeric(circular(bearing + rwrappedcauchy(n = 1,
                                                                                           mu = targeting_bias,
                                                                                           rho = concentration),
                                                                  modulo = "2pi"))
                                  location[1] <<- location + rgamma(1, shape = 1, scale = speed) * c(cos(bearing), sin(bearing)) #otherwise, move closer
                                }

                              } else { #if not targeting
                                bearing <<- as.numeric(circular(bearing, modulo = "2pi") + rwrappedcauchy(n = 1,
                                                                                                          mu = circular(turnBias, modulo = "2pi"),
                                                                                                          rho = concentration))
                                location[1] <<- location + rgamma(1, shape = 1, scale = speed) * c(cos(bearing), sin(bearing))
                              }

                              if (!is.na(bounds)){ #check for valid bounds
                                turnVarIncrease <- 0
                                while(! st_within(location, bounds, sparse = FALSE)[1,1]) {#check if out of bounds, if so . . .
                                  location[1] <<- st_point(path[[length(path)]][nrow(path[[length(path)]]),]) #reset location
                                  if(concentration-turnVarIncrease > 0.2) turnVarIncrease <- turnVarIncrease + 0.02 #relax directional concentration
                                  bearing <<- as.numeric(circular(bearing + rwrappedcauchy(n = 1,
                                                                                           mu = circular(0),
                                                                                           rho = concentration - turnVarIncrease),
                                                                  modulo = "2pi")) # get new bearing
                                  location[1] <<- location + speed * c(cos(bearing), sin(bearing)) #try moving again
                                }
                              }
                              path[[length(path)]] <<- st_linestring(rbind(path[[length(path)]], location[[1]]))
                              attr(path, "bbox") <<- setNames(unlist(bbox(as_Spatial(path))), c("xmin", "ymin", "xmax", "ymax"))
                              attr(location, "bbox") <<- setNames(unlist(bbox(as_Spatial(location))), c("xmin", "ymin", "xmax", "ymax"))
                            })
)

#' @title dForager Constructor
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
#' @description Constructor function for a Distance Discounting Biased Random Walk Forager class objects. Objects of this class have methods for finding targets, moving, and plotting. Movement switches between a biased random walk (with step lengths drawn from a gamma distribution and bearing drawn from a wrapped cauchy) and directed motion toward a target. Target selection is weighted toward closer targets.
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
#' @field concentration a numeric between 0 and one which gives the concentration of the wrapped cauchy distribution from which new bearings are drawn, where a 1 results in a turning angle equal to turnBias every step and a 0 results in a circular random uniform probability distribution
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
#' @title selector
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

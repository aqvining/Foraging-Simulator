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
                                                   patches = lapply(rep(2, times = 10), function(x) st_point(runif(x, min = -15, max = 15))) %>%
                                                     st_sfc %>% data.frame(geom = ., NAME = as.character(1:length(.)), MAX_VALUE = 10, VALUE = 8, REGEN = 1) %>%
                                                     st_sf() %>% st_buffer(1)
                                                   %>% st_set_crs(32610),
                                                   bounds = st_sfc(st_buffer(st_convex_hull(reduce(patches$geometry, c)),1), crs = st_crs(patches)), #default bounds created as convex hull around patches with a buffer of one
                                                   foragers = createForagers(3, bounds = bounds, speed = 2, quiet = TRUE) #default creation of forager within bounds
                             ){  "Set default values for variables that are not entered manually"
                               callSuper(..., foragers = foragers, patches = patches, bounds = bounds)
                             },
                             progress = function(){ "moves the entire environment forward one timestep by moving all foragers once and processing impact of foraging on patches"
                               for (forager in foragers[sample(length(foragers))]) {#operates on each forager in random order
                                 if (! forager$targeting & ! forager$foraging) forager$setTarget(patches) #If the forager does not currently have a target, look for one
                                 forager$move(bounds = bounds)
                                 if(forager$foraging) execute_forage(forager)
                               }
                               patches <<- patches %>% mutate(VALUE = logistic_growth(y = VALUE, max = MAX_VALUE, scale = REGEN))
                             },
                             execute_forage = function(forager) { #reduces the value of the forager's target patch, will eventually also increase energy of forager
                               updated_patch <- forager$target$NAME
                               extraction <- patches %>%            #Extracts 1/10th of VALUE from patch. Consider a more explicit model or variable to control the reduction
                                 st_drop_geometry() %>%
                                 filter(NAME == updated_patch) %>%
                                 select(VALUE) %>%
                                 unlist() * forager$efficiency           #extraction needs to be a numeric not a data.frame, but this feels hacky
                               patches <<- patches %>% mutate(VALUE = VALUE - (NAME %in% updated_patch) * extraction) #update value in patches
                               forager$energy <- forager$energy + extraction
                               forager$target <- patches %>% filter(NAME == updated_patch)
                               if(forager$giving_up_density > extraction) forager$foraging <- FALSE
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
                                      foragers[[1]]$path[[trials]] <<- st_cast(foragers[[1]]$location, "LINESTRING")[[1]]
                                      foragers[[1]]$bearing <<- as.numeric(rwrappedcauchy(n = 1, mu = circular(0), rho = 0))
                                      foragers[[1]]$visitSeq <<- rep("NA", times = foragers[[1]]$repeatAvoid)
                                    }
                                    callSuper()
                                    for (forager in foragers) {
                                      if (length(unique(forager$visitSeq[-c(1:forager$repeatAvoid)])) == nrow(patches) | nrow(forager$path[[1]]) >= 2500) { #end conditions for trial 1) Forager has been to all patches, 2) forager has made 2500 steps
                                        sequence <<- c(sequence, forager$visitSeq[-c(1:forager$repeatAvoid)])
                                        trials <<- trials + 1
                                      }
                                    }
                                  }
                                ))

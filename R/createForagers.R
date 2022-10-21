#' @title Forager Constructor Wrapper
#' @description A wrapper function to create a list of foragers from lists of parameter values. Handles single values if parameter is the same for all foragers. Assigns default values for empty fields.
#
#' @param numForagers A numeric giving the number of foragers to be created
#' @param type A character giving the type of forager to create. Currenly handles "Random" (default, parent class Forager) and "BRW" (brwForager)
#' @param bounds A simple features collection of polygons within which all foragers must be created
#' @param locations A list of simple features collections each containing a single st_POINT giving the initial location of a new forager. If NA or empty, foragers will be created at random locations within bounds
#' @param bearings A single numeric giving the initial bearing of all new foragers in radians or a list of numerics equal in length to the value of numForagers
#' @param speeds A single numeric giving the initial speed of all new foragers or a list of numerics equal in length to the value of numForagers
#' @param concentrations A single numeric giving the initial concentration of all new foragers or a list of numerics equal in length to the value of numForagers. Ignored unless type = "BRW"
#' @param persistences NUmeric vector with values between 0 and 1. Each element of this vector dictates the weight of an agents current direction relative to the direction of the agents current target when choosing a new direction
#' @param sights A single numeric giving the initial sight range of all new foragers or a list of numerics equal in length to the value of numForagers
#' @param giving_up_density the per-step energy return at which an agent will leave a foraging patch
#' @param efficiency the proportion of a patch that an agent extracts and adds to its energy each time step spend foraging
#' @param repeatAvoids A single numeric giving the number of unique patches a forager must visit before revisting a patch, or a list of numerics equal in length to the value of numForagers
#' @param quiet A logical indicating whether to suppress warnings about default values used to fill in empty fields
#' @param CRS A numeric giving the CRS code to assign to forager locations. If bounds argument is defined, the CRS of that object will be used and this argument will be ignored
#' @param turnBiases A single numeric giving the average turn angle of a forager in radians, or a list of numerics equal in length to the value of numForagers. Ignored unless type = "BRW"
#' @param choice_determinism A numeric vector with values between -1 and 1. For each agent, the corresponding value sets the randomness of its decisions, where -1 is nearly entirely random, 0 sets decision probabities directly proportional to the ratin of patch value to distance, and 1 nearly always selects the most valuable patch
#' @export
#' @importFrom circular circular rcircularuniform rwrappedcauchy
createForagers <- function(numForagers,
                           type = "Random",
                           bounds = NA,
                           locations = NA,
                           bearings = NA,
                           speeds = NA,
                           concentrations = NA,
                           persistences = NA,
                           sights = NA,
                           giving_up_density = NA,
                           efficiency = NA,
                           repeatAvoids = NA,
                           quiet = FALSE,
                           CRS = "NA",
                           turnBiases = NA,
                           choice_determinism = NA) {
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
  if (is.na(giving_up_density)) {
    if (! quiet) warning("No giving up density given, initial values set to 0.5")
    giving_up_density <- 0.5
  }
  if (is.na(efficiency)) {
    if (! quiet) warning("No efficiencies given, initial values set to 0.1")
    efficiency <- 0.1
  }
  if (is.na(repeatAvoids)) {
    if (! quiet) warning("No patch avoidance memory given. Initial values set to avoid 0 most recently visited patches")
    repeatAvoids <-  0
  }
  if (is.na(choice_determinism)){
    if (! quiet) warning("no choice determinism given. Initial values set to 0 (make choices proportional to attraction)")
    choice_determinism <- 0
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
  if (type == "BRW" | type == "dBRW"){
    if(is.na(concentrations)){
      if (! quiet) warning("No turning concentrations given, initial values set to 0.7")
      concentrations <- 0.7
    }
    if(is.na(persistences)){
      if (! quiet) warning("No directional persistences given, initial values set to 0.5")
      persistences <- 0.5
    }
    if(is.na(turnBiases)){
      if (! quiet) warning("No turning biases given, initial values set to 0")
      turnBiases <- 0
    }
  }
  if(!length(locations) == numForagers) stop("number of locations given must equal value of numForagers argument")
  locations <- lapply(locations, to_sf_point, crs = st_crs(bounds))
  foragers <- vector("list", length = numForagers)
  parameters <- st_sf(geom = reduce(locations,c), BEARING = bearings, SPEED = speeds, SIGHT = sights, GUD = giving_up_density, EFFICIENCY = efficiency, REPEATAVOID = repeatAvoids, DET = choice_determinism, concentration = concentrations, PERSISTENCE = persistences, BIAS = turnBiases)
  if(type == "Random") for (i in 1:numForagers) foragers[[i]] <- Forager(location=parameters$geom[i],
                                                                         bearing = parameters$BEARING[i],
                                                                         speed = parameters$SPEED[i],
                                                                         sight = parameters$SIGHT[i],
                                                                         giving_up_density = parameters$GUD[i],
                                                                         efficiency = parameters$EFFICIENCY[i],
                                                                         repeatAvoid = parameters$REPEATAVOID[i],
                                                                         choice_determinism = parameters$DET[i]
  )
  if(type == "BRW") for (i in 1:numForagers) foragers[[i]] <- BRWForager(location=parameters$geom[i],
                                                                         bearing = parameters$BEARING[i],
                                                                         speed = parameters$SPEED[i],
                                                                         sight = parameters$SIGHT[i],
                                                                         giving_up_density = parameters$GUD[i],
                                                                         efficiency = parameters$EFFICIENCY[i],
                                                                         repeatAvoid = parameters$REPEATAVOID[i],
                                                                         choice_determinism = parameters$DET[i],
                                                                         concentration = parameters$concentration[i],
                                                                         persistence = parameters$PERSISTENCE,
                                                                         turnBias = parameters$BIAS[i])
  if(type == "dBRW") for (i in 1:numForagers) foragers[[i]] <- DDBRWForager(location=parameters$geom[i],
                                                                            bearing = parameters$BEARING[i],
                                                                            speed = parameters$SPEED[i],
                                                                            sight = parameters$SIGHT[i],
                                                                            giving_up_density = parameters$GUD[i],
                                                                            efficiency = parameters$EFFICIENCY[i],
                                                                            repeatAvoid = parameters$REPEATAVOID[i],
                                                                            choice_determinism = parameters$DET[i],
                                                                            concentration = parameters$concentration[i],
                                                                            persistence = parameters$PERSISTENCE,
                                                                            turnBias = parameters$BIAS[i])
  return(foragers)
}

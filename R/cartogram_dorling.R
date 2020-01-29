#' @title Calculate Non-Overlapping Circles Cartogram
#' @description Construct a cartogram which represents each geographic region 
#' as non-overlapping circles (Dorling 1996).
#' @name cartogram_dorling
#' @param x SpatialPolygonsDataFrame, SpatialPointsDataFrame or an sf object
#' @param weight Name of the weighting variable in x
#' @param k Share of the bounding box of x filled by the larger circle
#' @param l Multiplying factor (in map units) to adjust the size of the circles.
#' Using l instead of k allows the cartograms to be made comparable with each other.
#' @param m_weight Circles' movements weights. An optional vector of numeric weights 
#' (0 to 1 inclusive) to 
#' apply to the distance each circle moves during pair-repulsion. A weight of 0 
#' prevents any movement. A weight of 1 gives the default movement distance. A 
#' single value can be supplied for uniform weights. A vector with length less 
#' than the number of circles will be silently extended by repeating the final 
#' value. Any values outside the range [0, 1] will be clamped to 0 or 1.
#' @param itermax Maximum iterations for the cartogram transformation. 
#' @return Non overlaping proportional circles of the same class as x.
#' @export
#' @references Dorling, D. (1996). Area Cartograms: Their Use and Creation. In Concepts and Techniques in Modern Geography (CATMOG), 59.
#' @examples
# Packages
#' library(sf)
#' library(rnaturalearth)
#' library(cartogram)
#' 
#' # Data inport and Handling
#' world <- ne_countries(returnclass = "sf")
#' crs <- "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
#' world <- st_transform(world, crs = crs)
#' world <- world[!is.na(world$pop_est),]
#' 
#' # Cartogram building
#' pop.dorling <- cartogram_dorling(world, "pop_est")
#' plot(st_geometry(pop.dorling), col ="red", border="white")
#' 
#' # Layout and Legend
#' title("World Population")

cartogram_dorling <- function(x, weight, k = 5, l = NULL, m_weight = 1, itermax= 1000) {
  UseMethod("cartogram_dorling")
}

#' @rdname cartogram_dorling
#' @importFrom sf st_is_longlat st_as_sf st_geometry st_coordinates st_geometry st_centroid st_crs
#' @importFrom packcircles circleRepelLayout
#' @export
cartogram_dorling.sf <- function(x, weight, k = 5, l = NULL, m_weight = 1, itermax= 1000){
  # proj or unproj
  if (sf::st_is_longlat(x)) {
    stop('Using an unprojected map. This function does not give correct centroids and distances for longitude/latitude data:\nUse "st_transform()" to transform coordinates to another projection.', call. = F)
  }
  # no 0 values
  x <- x[x[[weight]]>0,]
  # data prep
  dat.init <- data.frame(sf::st_coordinates(sf::st_centroid(sf::st_geometry(x))),
                         v = x[[weight]])
  surf <- (max(dat.init[,1]) - min(dat.init[,1])) *  (max(dat.init[,2]) - min(dat.init[,2]))
  # comparable or not
  if(!is.null(l) == TRUE){
    dat.init$v <- dat.init$v * l
  }
  if(!is.null(k) == TRUE){
    dat.init$v <- dat.init$v * (surf * k / 100) / max(dat.init$v)
  }
  # circles layout and radiuses
  res <- packcircles::circleRepelLayout(x = dat.init, xysizecols = 1:3,
                                        wrap = FALSE, sizetype = "area",
                                        maxiter = itermax, weights = m_weight)
  # sf object creation
  . <- sf::st_buffer(sf::st_as_sf(res$layout,
                                  coords =c('x', 'y'),
                                  crs = sf::st_crs(x)),
                     dist = res$layout$radius)
  sf::st_geometry(x) <- sf::st_geometry(.)
  return(x) 
}

#' @rdname cartogram_dorling
#' @export
cartogram_dorling.SpatialPolygonsDataFrame <- function(x, weight, k = 5, l = NULL, m_weight = 1, itermax = 1000){
  as(cartogram_dorling.sf(sf::st_as_sf(x), weight = weight, k = k, m_weight = m_weight, itermax = itermax), "Spatial")
}


#' @title Dorling Cartogram Legend
#' @description Build and display the legend of a Dorling cartogram
#' @name cartogram_dorling_legend
#' @param x dorling cartogram generetd with the function cartogram_dorling(). sf object.
#' @param var name of the vector containing the data 
#' @param pos position of the legend. a vector of two coordinates in map units (c(x, y)
#' @param col fill color.
#' @param border stroke color.
#' @param values.cex size of the values in the legend.
#' @param values.round rounding of numbers
#' @param lty line style
#' @param nb.circles number of circles in the legend
#' @param title.txt title of the legend.
#' @param title.cex size of the legend title.
#' @param title.font font of the legend title.
#' @export
#' @references Dorling, D. (1996). Area Cartograms: Their Use and Creation. In Concepts and Techniques in Modern Geography (CATMOG), 59.
#' @examples
# Packages
#' library(sf)
#' library(rnaturalearth)
#' library(cartogram)
#' 
#' # Data inport and Handling
#' world <- ne_countries(returnclass = "sf")
#' crs <- "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
#' world <- st_transform(world, crs = crs)
#' world <- world[!is.na(world$pop_est),]
#' 
#' # Cartogram building
#' pop.dorling <- cartogram_dorling(world, "pop_est")
#' plot(st_geometry(pop.dorling), col ="red", border="white")
#' 
#' # Layout and Legend
#' title("World Population")
#' cartogram_dorling_legend(x = pop.dorling, pos <- NULL, var = "pop_est", col = "red",
#'                          border = "black", lwd = 1, values.cex = 0.6, values.round = 0,
#'                          lty = 3, nb.circles =  4, title.txt = "Number of inhabitants.",
#'                         title.cex = 0.8, title.font = 2)
cartogram_dorling_legend <- function(x, var, pos = NULL, col = "white",
                                     border = "black",lwd = 1, values.cex = 0.6,
                                     values.round = 0, lty = 3, nb.circles = 4,
                                     title.txt = "Title of the legend",
                                     title.cex = 0.8,
                                     title.font = 2) {
  
  # Radii & Values
  v <- x
  st_geometry(v) <- NULL
  v <- v[,var]
  r <- sqrt(as.numeric(st_area(x))/pi)
  radii <- seq(from = max(r), to = min(r), length.out = nb.circles)
  sle <- radii * radii * pi
  values <- sle * max(v) / sle[1]
  
  # Positions
  par()$usr
  
  delta <- (par()$usr[2] - par()$usr[1]) / 50
  if(length(pos) != 2){
    pos <- c(par()$usr[1] + radii[1] + delta,par()$usr[3] + delta)
  }
  
  # Circles
  
  for(i in 1:nb.circles){
    # circles
    posx <- pos[1]
    posy <- pos[2] + radii[i]
    p <- st_sfc(st_point(c(posx,posy)))
    circle <- st_buffer(st_as_sf(p), dist = radii[i])
    plot(circle, col = col, border = border, lwd=lwd, add=T)
    # lines
    segments(posx, posy + radii[i], posx + radii[1] + radii[1]/10, col = border, lwd=lwd, lty = lty)
    # texts
    text(x = posx + radii[1] + radii[1]/5, y = posy + radii[i], 
         labels = formatC(round(values[i],values.round), big.mark = " ", format = "fg", digits = values.round), adj = c(0,0.5), cex = values.cex)
  }
  
  # Title
  text(x = posx - radii[1] ,y = posy + radii[1]*2 + radii[1]/3, title.txt,
       adj = c(0,0), cex = title.cex, font = title.font)
}
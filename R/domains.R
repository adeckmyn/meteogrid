compare.geodomain <- function(x, y){
### TRUE if they are equal, FALSE if they are not
### this is not exhaustive, but I am in a hurry!

### compare the projection: may have a different order!
  n1 <- names(x$projection)
  n2 <- names(y$projection)
  ndif <- setdiff(n1,n2)
  if (length(ndif) > 0) return(FALSE)
  if (!isTRUE(all.equal(x$projection, y$projection[n1]))) return(FALSE)

### now the main domain specifiers
### NOTE: clonlat OR SW-NE points may be known
  checklist <- c("nx", "ny", "dx", "dy", "clonlat")
  if (any(is.null(x[checklist]))) x <- DomainExtent(x)
  if (any(is.null(y[checklist]))) y <- DomainExtent(y)
  if (!isTRUE(all.equal(x[checklist], y[checklist]))) return(FALSE)
  else return(TRUE)
  # $ex and $ey may be different, so ignore them
}

###########################
### CREATE NEW DOMAINS  ###
###########################

Make.domain <- function(projtype="lambert", clonlat, nxny, dxdy, exey=NULL,
                        reflat=clonlat[2], reflon=clonlat[1], tilt=0,
                        earth=list(R=6371229)){
  if (length(dxdy)==1) dxdy <- rep(dxdy,2)
  if (length(nxny)==1) nxny <- rep(nxny,2)
  if (!is.null(exey) && length(exey)==1) exey <- rep(exey,2)

  if (projtype %in% c("lcc", "lambert")) {
### Lambert (as used in ALADIN: only 1 reference latitude)
    projection <- list(proj="lcc", lon_0=reflon, lat_1=reflat, lat_2=reflat)
  } else if (projtype %in% c("merc", "mercator")){
### Mercator
    projection <- list(proj="merc",lon_0=clonlat[1])
  } else if (projtype %in% c("omerc", "tmerc")){
### Rotated & tilted Mercator
### as used in ALADIN: reflon=clon, reflat=clat !!!
### we simply *ignore* reflon and reflat (coming from FA file they are 0 !)
#    if (reflat != clonlat[2] || reflon != clonlat[1]) warning('This domain is not ALADIN-compatible!')
    rlon <- clonlat[1]
    rlat <- clonlat[2]
    if (abs(rlat) < 0.01 && abs(tilt) < 0.01) {
# not necessary: this is OK now in "omerc"
      projection <- list(proj="merc",lon_0=rlon)
    } else if (abs(abs(tilt)-90) < 1.0E-5) {
# this would crash in "omerc" (alpha ~ 0)
      projection <- c(proj = "tmerc", lon_0 = rlon, lat_0 = rlat)
    } else if (abs(tilt) < 1.E-5) {
      projection <- c(proj = "somerc", lon_0 = rlon, lat_0 = rlat)
    } else if (tilt>0) {
      projection <- list(proj = "omerc", lonc = rlon,
                            lat_0 = rlat, alpha = -90 + tilt, no_rot=NA)
    } else {
      projection <- list(proj = "omerc", lonc = rlon,
                            lat_0 = rlat, alpha = 90 + tilt, no_rot=NA)
    }
  } else if (projtype %in% c("latlong", "lalo", "lonlat")) {
    projection=list(proj="latlong")
  } else if (projtype %in% c("ob_tran", "RotLatLon")){
    projection <- list(proj="ob_tran","o_proj"="latlong",
                       "o_lat_p"=-reflat,"o_lon_p"=0,"lon_0"=reflon)
  } else if (projtype %in% c("stere", "stereographic")) {
    projection <- list(proj="stere",lon_0=reflon,lat_0=reflat)
  } else {
    stop("Unknown projection.")
  }
  projection <- c(projection, earth)
### project the centre point

#  cxy <- project(list(x=clonlat[1], y=clonlat[2]), proj=projection)

#  SW0 <- c(cxy$x,cxy$y) - dxdy*(nxny - 1)/2
#  NE0 <- SW0 + dxdy*(nxny - 1)

#  lims <- project(list(x=c(SW0[1],NE0[1]),y=c(SW0[2],NE0[2])), proj=projection, inv=TRUE)
#  SW <- c(lims$x[1],lims$y[1])
#  NE <- c(lims$x[2],lims$y[2])
### and the output is...
  result <- list(projection=projection, nx=nxny[1], ny=nxny[2], dx=dxdy[1], dy=dxdy[2],
                 clonlat=clonlat)
  if (!is.null(exey)) result <- c(result, ex=exey[1], ey=exey[2])
  class(result) <- "geodomain"
  result
}

Make.domain.RLL <- function(Lon1, Lat1, SPlon, SPlat, SPangle=0, nxny, dxdy){
### This is for Rotated LatLon as used by Hirlam: central meridian is vertical.
### In the future, this should be merged with Make.domain
### but that is not trivial: here, the centre point is of no real consequence
### and you have to define the rotated South Pole
### RLL projection uses radians, so we have to rescale dxdy at the end
### BUT that is rather ugly...
  if (length(dxdy)==1) dxdy <- rep(dxdy, 2)
  if (length(nxny)==1) nxny <- rep(nxny, 2)
  cxy <- (c(Lon1, Lat1) + (nxny -1) * dxdy / 2) / 180 * pi
  projection <- list(proj="ob_tran", "o_proj"="latlong",
                       "o_lat_p"=-SPlat, "o_lon_p"=0, "lon_0"=SPlon)
  cll <- project(cxy[1], cxy[2], proj=projection, inv=TRUE)
  
#  Lon2 <- Lon1 + (nxny[1]-1)*dxdy[1]
#  Lat2 <- Lat1 + (nxny[2]-1)*dxdy[2]
#  RR <- project(list(x=c(Lon1, Lon2)*pi/180, y=c(Lat1, Lat2)*pi/180), proj=projection, inv=TRUE)
#  SW <- c(RR$x[1], RR$y[1])
#  NE <- c(RR$x[2], RR$y[2])
  
  result <- list(projection=projection, nx=nxny[1], ny=nxny[2],
		 clonlat=c(cll$x, cll$y), dx=dxdy[1]/180*pi, dy=dxdy[2]/180*pi)
  class(result) <- "geodomain"
  result
}

#################################
### Domain points and borders ###
#################################

DomainExtent <- function(geo) {
### We look for the extreme LatLon values of a domain domain
### these are useful to select the right part from the world map
### first find the SW and NE point, then draw a rectangular box
### and project it back to LatLon
### You could use DomainPoints in stead, but then you'd be projecting
### all domain points, which is a bit of an overkill...
### old domain definitions may have SW and NE points, but newer ones have central point only
  domain <- as.geodomain(geo)

  if (!is.null(domain$clonlat) && !is.null(domain$dx) && !is.null(domain$dy)) {
    dx <- domain$dx
    dy <- domain$dy
    clonlat <- domain$clonlat
    cxy <- project(domain$clonlat, proj=domain$projection)
    x0  <- cxy$x - dx*(domain$nx -1)/2
    x1  <- cxy$x + dx*(domain$nx -1)/2
    y0  <- cxy$y - dy*(domain$ny -1)/2
    y1  <- cxy$y + dy*(domain$ny -1)/2
  } else {
    xy <- project(list(x=c(domain$SW[1], domain$NE[1]),y=c(domain$SW[2], domain$NE[2])),
                proj=domain$projection)
    x0 <- xy$x[1]
    y0 <- xy$y[1]
    x1 <- xy$x[2]
    y1 <- xy$y[2]
    dx <- (x1-x0)/(domain$nx-1)
    dy <- (y1-y0)/(domain$ny-1)
    xc <- (x0+x1)/2
    yc <- (y0+y1)/2
    cll <- project(list(x=xc,y=yc), proj=domain$projection, inv=TRUE)
    clonlat <- c(cll$x, cll$y)
  }
  border_bottom <- project(list(x=seq(x0, x1, length.out=domain$nx), y=rep(y0, domain$nx)),
                      proj=domain$projection, inv=TRUE)
  border_right <- project(list(x=rep(x1, domain$ny), y=seq(y0, y1, length.out=domain$ny)),
                      proj=domain$projection, inv=TRUE)
  border_top <- project(list(x=seq(x0, x1, length.out=domain$nx), y=rep(y1, domain$nx)),
                      proj=domain$projection, inv=TRUE)
  border_left <- project(list(x=rep(x0, domain$ny), y=seq(y0, y1, length.out=domain$ny)),
                      proj=domain$projection, inv=TRUE)

  borders <- project(list(x=c(seq(x0, x1, length.out=domain$nx), rep(x1,domain$ny),
                              seq(x0, x1, length.out=domain$nx), rep(x0,domain$ny)),
                          y=c(rep(y0, domain$nx), seq(y0,y1, length.out=domain$ny),
                              rep(y1, domain$nx), seq(y0,y1, length.out=domain$ny))),
                      proj=domain$projection, inv=TRUE)

  ### FIX ME: should we make sure that, before checking, all borders are in [-180, 180] ?
  ### look for *large* negative jumps (e.g. 360 -> 0)

  t1 <- which(diff(border_top$x) < -300)
  if (length(t1) == 1) { # there is a date line jump
    tail(border_top$x, -t1) <- tail(border_top$x, -t1) + 360
    if (any(border_top$x > 360)) border_top$x <- border_top$x - 360
  } else if (length(t1) > 1) {
    stop("Multiple longitude jumps detected. Don't know how to handle.")
  }

  t2 <- which(diff(border_bottom$x) < -300)
  if (length(t2) == 1) { # there is a date line jump
    tail(border_bottom$x, -t2) <- tail(border_bottom$x, -t2) + 360
    if (any(border_bottom$x > 360)) border_bottom$x <- border_bottom$x - 360
  } else if (length(t2) > 1) {
    stop("Multiple longitude jumps detected. Don't know how to handle.")
  }

  lonlim <- range(c(border_top$x, border_right$x, border_bottom$x, border_left$x))
  latlim <- range(c(border_top$y, border_right$y, border_bottom$y, border_left$y))
  if (lonlim[2] > 180) wrap <- c(0,360)
  else wrap <- FALSE

  list(lonlim = lonlim , latlim = latlim,
       clonlat=clonlat,
       x0=x0, y0=y0, x1=x1, y1=y1, dx=dx, dy=dy, nx=domain$nx, ny=domain$ny,
       NE=c(border_top$x[domain$nx], border_top$y[domain$nx]),
       SW=c(border_bottom$x[1], border_bottom$y[1]),
       wrap = wrap)
}

##########################################################
DomainPoints <- function (geo, type="lalo"){
### return lat's and lon's of all domain points (or leave in projection if type = "xy")
  domain <- as.geodomain(geo)
 
  if (!is.null(domain$clonlat) && !is.null(domain$dx) && !is.null(domain$dy)) {
    cxy <- project(domain$clonlat, proj=domain$projection)
    xy <- data.frame(x= cxy$x + c(-1, +1) * domain$dx * (domain$nx-1)/2 ,
                     y= cxy$y + c(-1, +1) * domain$dy * (domain$ny-1)/2) 
  } else {
    lalo <- list(x=c(domain$SW[1],domain$NE[1]),y=c(domain$SW[2],domain$NE[2]))
    xy <- project(lalo, proj = domain$projection)
  }
  xydomain <- expand.grid(x = seq(xy$x[1], xy$x[2], length.out = domain$nx),
                          y = seq(xy$y[1], xy$y[2], length.out = domain$ny))
  if (type == "lalo") {
    lalolist <- project(xydomain, proj = domain$projection, inv = TRUE)
    list(lon = matrix(lalolist$x, ncol = domain$ny, nrow = domain$nx),
         lat = matrix(lalolist$y, ncol = domain$ny, nrow = domain$nx))
  }
  else if (type == "xy") {
    list(x=matrix(xydomain$x, ncol = domain$ny, nrow = domain$nx),
         y=matrix(xydomain$y, ncol = domain$ny, nrow = domain$nx))
  } else {
    print("Unknown type.")
  }
}

DomainCorners <- function(geo) {
  domain <- as.geodomain(geo)

  if (!is.null(domain$clonlat) && !is.null(domain$dx) && !is.null(domain$dy)) {
    cxy <- project(domain$clonlat, proj=domain$projection)
    x0  <- cxy$x - domain$dx*(domain$nx -1)/2
    x1  <- cxy$x + domain$dx*(domain$nx -1)/2
    y0  <- cxy$y - domain$dy*(domain$ny -1)/2
    y1  <- cxy$y + domain$dy*(domain$ny -1)/2
    SW <- project(x0, y0, proj=domain$projection, inv=TRUE)
    NE <- project(x1, y1, proj=domain$projection, inv=TRUE)
  } else {
    xy <- project(list(x=c(domain$SW[1], domain$NE[1]),
                       y=c(domain$SW[2], domain$NE[2])),
                  proj=domain$projection)
    x0 <- xy$x[1]
    y0 <- xy$y[1]
    x1 <- xy$x[2]
    y1 <- xy$y[2]
    SW <- domain$SW
    NE <- domain$NE
  }
  SE <- project(x1, y0, proj=domain$projection, inv=TRUE)
  NW <- project(x0, y1, proj=domain$projection, inv=TRUE)
  list(SW=SW, NE=NE, SE=SE, NW=NW)
}

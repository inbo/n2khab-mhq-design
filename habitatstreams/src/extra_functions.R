#' @title split lines
#' @description Splits lines longer than a given threshold into the minimum number of pieces to all be under the given threshold.
#' @param lines data.frame of class sf with LINESTRING sfc column.
#' @param max_length maximum segment length to return
#' @param id name of ID column in data.frame
#' @return only the split lines.
#' @importFrom dplyr group_by ungroup filter left_join select rename mutate
#' @export
#'
split_lines <- function(input_lines, max_length, id = "ID") {
  geom_column <- attr(input_lines, "sf_column")
  
  input_crs <- sf::st_crs(input_lines)
  
  input_lines[["geom_len"]] <- sf::st_length(input_lines[[geom_column]])
  
  attr(input_lines[["geom_len"]], "units") <- NULL
  input_lines[["geom_len"]] <- as.numeric(input_lines[["geom_len"]])
  
  too_long <- filter(select(input_lines, id, geom_column, geom_len), geom_len >= max_length)
  
  rm(input_lines) # just to control memory usage in case this is big.
  
  too_long <- mutate(too_long,
                     pieces = ceiling(geom_len / max_length),
                     piece_len = (geom_len / pieces),
                     fID = 1:nrow(too_long))
  
  split_points <- sf::st_set_geometry(too_long, NULL)[rep(seq_len(nrow(too_long)), too_long[["pieces"]]),]
  
  split_points <- mutate(split_points, split_fID = row.names(split_points)) %>%
    select(-geom_len, -pieces) %>%
    group_by(fID) %>%
    mutate(ideal_len = cumsum(piece_len)) %>%
    ungroup()
  
  coords <- data.frame(sf::st_coordinates(too_long[[geom_column]]))
  rm(too_long)
  
  coords <- rename(coords, fID = L1) %>% mutate(nID = 1:nrow(coords))
  
  split_nodes <- group_by(coords, fID) %>%
    # First calculate cumulative length by feature.
    mutate(len  = sqrt(((X - (lag(X)))^2) + (((Y - (lag(Y)))^2)))) %>%
    mutate(len = ifelse(is.na(len), 0, len)) %>%
    mutate(len = cumsum(len)) %>%
    # Now join nodes to split points -- this generates all combinations.
    left_join(select(split_points, fID, ideal_len, split_fID), by = "fID") %>%
    # Calculate the difference between node-wise distance and split-point distance.
    mutate(diff_len = abs(len - ideal_len)) %>%
    # regroup by the new split features.
    group_by(split_fID) %>%
    # filter out na then grab the min distance
    filter(!is.na(diff_len) & diff_len == min(diff_len)) %>%
    ungroup() %>%
    # Grab the start node for each geometry -- the end node of the geometry before it.
    mutate(start_nID = lag(nID),
           # need to move the start node one for new features.
           new_feature = fID - lag(fID, default = -1),
           start_nID = ifelse(new_feature == 1, start_nID + 1, start_nID)) %>%
    # Clean up the mess
    select(fID, split_fID, start_nID, stop_nID = nID, -diff_len, -ideal_len, -len, -X, -Y)
  
  split_nodes$start_nID[1] <- 1
  
  split_points <- left_join(split_points, select(split_nodes, split_fID, start_nID, stop_nID), by = "split_fID")
  
  new_line <- function(start_stop, coords) {
    sf::st_linestring(as.matrix(coords[start_stop[1]:start_stop[2], c("X", "Y")]))
  }
  
  split_lines <- apply(as.matrix(split_points[c("start_nID", "stop_nID")]),
                       MARGIN = 1, FUN = new_line, coords = coords)
  
  split_lines <- st_sf(split_points[c(id, "split_fID")], geometry = st_sfc(split_lines, crs = input_crs))
  
  return(split_lines)
}


#' @title Split Lines
#' @description to be filled
#' @param spobj sp object
#' @param dist distance in units of coords
#' @param start start from the start (T) or the end (F)
#' @return spatialLinesDataFrame
#' @import sp
#' @import plyr
#' @export
splitLines<- function(spobj, dist, start = T, sf = F){
  xydf<-coordBuild(spobj)
  if (start == F){
    xydf<-xydf[rev(rownames(xydf)),]
  }
  spoints <- split(xydf, dist)
  linelist <- list()
  lineslist <- list()
  id <- 1
  if(!sf) {
    j <- 1
    for(i in 1:(nrow(spoints)-1)){
      linelist[j] <- Line(spoints[c(i, i + 1), c(1:2)])
      j = j + 1
      if(spoints[i+1,3] == 1){ 
        lineslist[id]<-Lines(linelist, ID = id)
        id = id+1
        linelist<-list()
        j = 1
      }
    }
    return(SpatialLinesDataFrame(SpatialLines(lineslist), data = data.frame(id = 0:(length(lineslist)-1))))
  } else {
    start <- 1
    for(i in 1:(nrow(spoints)-1)){
      if(spoints[i+1,3] == 1){ 
        lineslist[[id]] <- sf::st_linestring(as.matrix(spoints[c(start:(i + 1)), c(1:2)], ncol = 2))
        id <- id + 1
        start <- i + 1
      }
    }
    return(sf::st_sf(id = 1:length(lineslist), geom = sf::st_sfc(lineslist)))
  }
}

#' @title build x/y coordinates data.frame
#' @description to be filled
#' @param spobj sp or sf object
#' @return x/y data.frame
#' @import sp
#' @import plyr
coordBuild <- function(spobj){
  if("LINESTRING" %in% class(spobj)) {
    return(setNames(data.frame(sf::st_zm(sf::st_coordinates(spobj))), c("x", "y")))
  }
  if(class(spobj) %in% c("SpatialLinesDataFrame",    "SpatialLines")){
    coords <- lapply(spobj@lines, function(x) lapply(x@Lines, function(y) y@coords))
    coords <- ldply(coords, data.frame)
    names(coords) <- c("x","y")
    return(coords)
  }
  if(class(spobj) %in% c("SpatialPolygonsDataFrame", "SpatialPolygons")){
    coords <- lapply(spobj@polygons, function(x) lapply(x@Polygons, function(y) y@coords))
    coords <- ldply(coords, data.frame)
    names(coords) <- c("x","y")
    return(coords)
  }
  if(class(spobj) == "data.frame"){
    if(all(c("x", "y") %in% tolower(names(spobj)))){
      return(spobj[c("x", "y")])
    }
    else{
      stop("Dataframe provided does not have x, y columns")
    }
  }
  stop("Class of spatial argument is not supported. Need SpatialLinesDataFrame or SpatialPolygonsDataFrame")
}

#' @title reads in xy data frame of coordinates and a distance value to split line on
#' @description to be filled
#' @param xydf dataframe with x and y coords
#' @param dist distance in units of coords
#' @return split version of xydf
#' @import sp
#' @import plyr
#' @importFrom utils tail
#'
split <- function(xydf, dist){
  modck <- function(change, mod){
    if(change<0){
      return(mod*-1)
    }
    else{
      return(mod)
    }
  }
  x <- c()
  y <- c()
  end <- c()
  rem <- 0
  for(i in 1:nrow(xydf)){
    if(i == 1){
      x <- c(x, xydf$x[i])
      y <- c(y, xydf$y[i])
      end <- c(end,1)
    }
    if(i != 1){
      cx <- xydf$x[i] - xydf$x[i-1]
      cy <- xydf$y[i] - xydf$y[i-1]
      if(cx & cy != 0){
        len <- sqrt((cx^2) + (cy^2)) + rem
        segs <- len %/% dist
        if(segs == 0){
          rem <- len
        }
        else{
          m <- cx/cy
          ymod <- dist / (sqrt((m^2)+1))
          xmod <- ymod * abs(m)
          yremsub <- rem / (sqrt((m^2)+1))
          xremsub <- yremsub * abs(m)
          xmod <- modck(cx, xmod)
          ymod <- modck(cy, ymod)
          xremsub <- modck(cx, xremsub)
          yremsub <- modck(cy, yremsub)
          xnew <- seq(xydf$x[i-1] - xremsub, xydf$x[i-1] + (xmod * segs), by = xmod)[-1]
          ynew <- seq(xydf$y[i-1] - yremsub, xydf$y[i-1] + (ymod * segs), by = ymod)[-1]
          if(length(xnew) != length(ynew)){
            if(abs(length(xnew) - length(ynew)) > 1) stop("Error found in new sequence. Code needs to be reviewed...")
            if(length(xnew) < length(ynew)){
              xnew <- c(xnew, xydf$x[i-1] + (xmod * segs))
            } else {
              ynew <- c(ynew, xydf$y[i-1] + (ymod * segs))
            }
          }
          rem<-sqrt((xydf$x[i] - tail(xnew,1))^2 + (xydf$y[i] - tail(ynew,1))^2)
          x <- c(x, xnew)
          y <- c(y, ynew)
          end <- c(end, rep(1, length(xnew)))
        }
      }
      if(cx != 0 & cy == 0){
        len <- cx + rem
        segs <- len %/% dist
        if(segs == 0){
          rem <- len
        }
        else{
          xmod <- dist
          ymod <- 0
          xmod <- modck(cx, xmod)
          xremsub <- modck(cx, xremsub)
          yremsub <- 0     
          xnew <- seq(xydf$x[i-1] - rem, xydf$x[i-1] + (xmod * segs), by = xmod)[-1]
          ynew <- rep(xydf$y[i-1], segs)
          rem <- xydf$x[i] - tail(xnew,1)
          x <- c(x, xnew)
          y <- c(y, ynew)
          end <- c(end, rep(1, length(xnew)))
        }
      }
      if(cx == 0 & cy != 0){
        len <- cy + rem
        segs <- len %/% dist
        if(segs == 0){
          rem <- len
        }
        else{
          xmod <- 0
          ymod <- dist
          xmod <- modck(cx, xmod)
          xremsub <- modck(cx, xremsub)
          yremsub <- 0    
          xnew <- rep(xydf$x[i-1], segs) 
          ynew <- seq(xydf$y[i-1] - rem, xydf$y[i-1] + (ymod * segs), by = ymod)[-1]
          rem <- xydf$y[i] - tail(ynew,1)
          x <- c(x, xnew)
          y <- c(y, ynew)
          end <- c(end, rep(1, length(ynew)))
        }
      }
      x <- c(x, xydf$x[i])
      y <- c(y, xydf$y[i])
      if(i != nrow(xydf)){
        end <- c(end, 0)    
      }
      else{
        end <- c(end, 1) 
      }
    }
  }
  return(data.frame(x = x, y = y, end = end))
}


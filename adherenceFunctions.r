#library(ntlpGlobals)
library(xts)

## lenOfTherapy assumes valid data, returns integer
lenOfTherapy <- function(pmc){
    last(pmc[,"fillDay"]) - first(pmc[,"fillDay"]) + last(pmc[,"days"])
}

## lenOfTherapy(pmc)

## persistence returns T/F, T if lenOfTherapy > threshold, F otherwise
## or a % relative to the threshold.  If threshold is not supplied return NA
## would it be usefull to have a % of threshold value too?
persistence <- function(pmc, threshold=NA, percent=FALSE){
    if(is.na(threshold) | !is.numeric(threshold)) return(NA)
    if(percent){
        lenOfTherapy(pmc) / threshold
    } else {
        lenOfTherapy(pmc) >= threshold
    }
}

## persistence(pmc,threshold=2000,percent=TRUE)
## persistence(pmc,threshold=2000,percent=F)
## persistence(pmc,1000,TRUE)
## persistence(pmc,1000,F)

## Days of coverage typically is the total days supply of med's the patient had on hand
## during a year (or alternate interval).  This is complicated somewhat by not having
## dates (and only number of days from index date) but not much.
## There are a number of ways this can be calculated:  worst case coverage over entire period,
## average case, best case, etc.
## Also can return a integer (as typically seen) or a percentage over the period.
## Also, may want to drop the last period since we don't know the end day.
daysOfCover <- function(pmc,termLen=1.0){
    ## maybe do this by year only, for each year averages?
    p1 <- pmc
    dayNum <- pmc[,"fillDay"] - first(pmc[,"fillDay"])
    ## cumDays <- cumsum(pmc[,"days"])
    ## cumDays <- c(NA,cumDays[-length(cumDays)])
    p1 <- cbind(pmc,dayNum)
    endDay <- p1[dim(p1)[1],"dayNum"]
    stDay <- endDay - ceiling(365*termLen)
    sel <- p1[,"dayNum"] > stDay
    p1 <- p1[sel,]
    p1[,"dayNum"] <- p1[,"fillDay"] - first(p1[,"fillDay"])
    cumDays <- cumsum(p1[,"days"])
    cumDays <- c(NA,cumDays[-length(cumDays)])
    p1 <- cbind(p1,cumDays)
    last(p1[,"cumDays"])
}

## daysOfCover(pmc)
## daysOfCover(pmc,termLen=2.0)

## medGaps computes the gap from inferred depletion date to refill date.
## Doesn't address changes in dossage.
medGaps <- function(pmc,evalType=c("value","count"),
                    evalCase=c("sum","mean","max","min","median"),
                    byYear=TRUE){
    totalDays <- last(pmc[,"fillDay"]) - first(pmc[,"fillDay"])
    gaps <- NULL
    p1 <- pmc
    repeat{
        dayNum <- p1[,"fillDay"] - first(p1[,"fillDay"])
        cumDays <- cumsum(p1[,"days"])
        cumDays <- c(NA,cumDays[-length(cumDays)])
        p1 <- cbind(p1,dayNum,cumDays)
        i1 <- which(!(p1[,"dayNum"] <= p1[,"cumDays"]))[1] # find first gap
        if(is.na(i1)) break ## no gaps then done
        g1 <- p1[i1,"dayNum"] - p1[i1,"cumDays"]
        gaps <- c(gaps,g1)
        p1 <- p1[i1:dim(p1)[1],]
        p1[,"dayNum"] <- p1[,"dayNum"] - p1[1,"dayNum"]
        p1 <- p1[,1:4]
    }
    if(evalType[1]=="count"){  # use number of gaps
        if(byYear){
            length(gaps)/totalDays*365
        } else {
            length(gaps)
        }
    } else { # use values
        if(byYear){
            do.call(evalCase[1],list(gaps))/totalDays*365
        } else {
            do.call(evalCase[1],list(gaps))
        }
    }
}

## medGaps(pmc,evalCase="min",byYear=F)
## medGaps(pmc,evalCase="mean",byYear=F)
## medGaps(pmc,evalCase="sum")
## medGaps(pmc,evalType="count")
## medGaps(pmc,evalType="count",byYear=F)
## medGaps(pmc) # default is to sum gap days normalize by year

## medPosRatio computes MPR, medication possession ratio given pmc.
## types: lot: length of therapy, or "last" for the last
## lastTerm*365 days on record.
medPosRat <- function(pmc,type=c("lot","last"),termLen=1.0){
    p1 <- pmc
    dayNum <- p1[,"fillDay"] - first(p1[,"fillDay"])
    p1 <- cbind(p1,dayNum)
    if(type[1] == "last"){ # crop data to a amount of year(s)
        endDay <- p1[dim(p1)[1],"dayNum"]
        stDay <- endDay - ceiling(365*termLen)
        stDay
        sel <- p1[,"dayNum"] > stDay
        p1 <- p1[sel,]
        p1[,"dayNum"] <- p1[,"fillDay"] - first(p1[,"fillDay"])
    }
    cumDays <- cumsum(p1[,"days"])
    cumDays <- c(NA,cumDays[-length(cumDays)])
    p1 <- cbind(p1,cumDays)
    p1
    p1[dim(p1)[1],"cumDays"]/p1[dim(p1)[1],"dayNum"]
}

## medPosRat(pmc)
## medPosRat(pmc,type="last")
## medPosRat(pmc,type="last",lastTerm=1.75)

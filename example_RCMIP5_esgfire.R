
# library(RCMIP5)

# need to do source('sourceall.R') for now


path <- "G:/Yannick"
RANGE1 <- 1991:2010
RANGE2 <- 2081:2100

varlist <- c("tas")
model <- "HadGEM2-CC"
experiments <- c("historical","rcp45", "rcp85")
tf <- tempfile()

first=TRUE
for(v in varlist) {
    print("================================================")
    print(v)
    print("================================================")
    for(ex in experiments) {
        print("---------------")
        print(ex)
        print("---------------")
        d <- loadEnsemble(variable=v, model=model, experiment=ex, path=path, ensemble="r1i1p1", verbose=T)
        print(summary(d))
        if (is.null(d)) next
        d <- filterDimensions(d, years=c(RANGE1, RANGE2), verbose=T)
        d <- makeMonthlyStat(d, verbose=T)
        print(summary(d))
        
        df <- as.data.frame(d)
        print(summary(df))
        write.table(df, tf, row.names=F, sep=",", col.names=first, append=!first)
        first=FALSE
    }
}

results <- read.csv(tf)

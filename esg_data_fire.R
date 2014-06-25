# Script to process ESG files into CSV annual summary format
# B.B.-L. and C.H. for Yannick fire project
# March 2014
# combines both land and ocean data into one script with latitude and longitude bands

# Important variable definitions, esp. data source & destination
SCRIPTNAME		<- "esg_data_fire.R"
INPUT_DIR		<- "G:/Yannick"
OUTPUT_DIR		<- "G:/Yannick/output/"
LOG_DIR			<- "logs/"
SEPARATOR		<- "-------------------"

# To use:
#   1. Set input and output directories, above
#   2. Set year ranges, below
#   3. Set scenario, below
#	4. Call process_variable

YEAR_RANGE1			<- 2001:2005
YEAR_RANGE2			<- 2096:2100	# these should be same length
SCENARIO            <- "rcp85"      # plus any historical files will also be processed
VARIABLES           <- c( "tas")
MODELS              <- c( "CMCC-CESM", "CanESM2")

DATAFREQ_ANNUAL     <- "annual"
DATAFREQ_MONTHLY    <- "monthly"
MONTHS_PER_YEAR   	<- 12

LAT_NAME      		<- "lat"
LAT_NAME2       	<- "j"
LAT_NAME3         <- "rlat"
LAT_NAME4         <- "x"
LEVEL_DEPTH     	<- 100 #20000 # from 10000000 to 20000 the upper limits of the troposphere.
LEVEL_NAME      	<- "lev" # change to lev for ocean data or plev for atmo data
LON_NAME      		<- "lon"
LON_NAME2       	<- "i"
LON_NAME3         <- "rlon"
LON_NAME4         <- "y"

INVFILE				<- paste0( OUTPUT_DIR, "inventory.csv" )
TIME_NAME       	<- "time"

# -----------------------------------------------------------------------------
# Load requested libraries
loadlibs <- function( liblist ) {
	printlog( "Loading libraries..." )
	loadedlibs <- vector()
	for( lib in liblist ) {
		printlog( "Loading", lib )
		loadedlibs[ lib ] <- require( lib, character.only=T )
		if( !loadedlibs[ lib ] )
			warning( "this package is not installed!" )
	}
	invisible( loadedlibs )
} # loadlibs

# -----------------------------------------------------------------------------
# Print dimensions of data frame
printdims <- function( d, dname=deparse( substitute( d ) ) ) {
	stopifnot( is.data.frame( d ) )
	printlog( dname, "rows =", nrow( d ), "cols =", ncol( d ) )
} # printdims

# -----------------------------------------------------------------------------
# Time-stamped output function
printlog <- function( msg, ..., ts=TRUE, cr=TRUE ) {
	if( ts ) cat( date(), " " )
		cat( msg, ... )
	if( cr ) cat( "\n")
}

# -----------------------------------------------------------------------------
# We keep an inventory of processed files
invfile <- function( fname="", model="", scenario="", ensemble="", units="", status="OK", newfile=FALSE ) {
	if( newfile | !file.exists( INVFILE ) ) {
		cat( "Date,File,Model,Years,Scenario,Ensemble,Units,Status\n", file=INVFILE )
	}
	if( fname!= "" | status != "") {
		cat( date(), fname, model, paste( range( YEAR_RANGE1 ), collapse="-" ), scenario, 
		ensemble, units, status, "\n", file=INVFILE, sep=",", append=T )
	}
}

# -----------------------------------------------------------------------------
# Get names of dimensions for a particular variable
read_dimension_names <- function( ncid, variable ) {
	dnames <- c()
	for( i in 1:ncid$nvars ) {
		if( ncid$var[[ i ]]$name == variable ) {
			for( j in 1:ncid$var[[ i ]]$ndims ) {
				dnames[ j ] <- ncid$var[[ i ]]$dim[[ j ]]$name
			}
			return( dnames )  
		}
	} # for
	stop( "*** Variable", variable, "not found!!" )
} # read_dimension_names

# -----------------------------------------------------------------------------
# Figure out whether there are levels in the file, and if so, how many we need to use
# to approximately reach a particular depth
determine_levels <- function( ncid, dnames, depth ) {
	levels_to_read <- 0
	if( LEVEL_NAME %in% dnames ) {
		levels_to_read <- which.min( abs( get.var.ncdf( ncid, LEVEL_NAME )-depth ) )
	}
	printlog( "Levels to read for", LEVEL_DEPTH, "=", levels_to_read )
	return( levels_to_read )
} # determine_levels

# -----------------------------------------------------------------------------
# Read a single time point, averaging over levels if necessary
read_timepoint <- function( ncid, variable, levels_to_avg, startdata, countdata, levindex ) {
	# printlog( "read_timepoint, startdata=", startdata, ", countdata=", countdata, ", variable=", variable, "levindex=", levindex )
	if( levels_to_avg >0 ) {
		dsum <- data.frame()
		for( i in 1:levels_to_avg ) {
			# printlog( "--     reading level", i )
			startdata[ levindex ] <- i
			dtemp <- get.var.ncdf( ncid, variable, start=startdata, count=countdata )
			if( i==1 )
			dsum <- dtemp
			else
			dsum <- dsum + dtemp 
		}
		d <- dsum / levels_to_avg
	} else {
		d <- get.var.ncdf( ncid, variable, start=startdata, count=countdata ) 
	}

	return( d )
} # read_timepoint

# -----------------------------------------------------------------------------
# The workhorse: read a file's data, fill data frame, write out
process_data <- function( fn, variable, beginyear, endyear, datafreq, year_range ) {

	stopifnot( datafreq==DATAFREQ_MONTHLY )

	results <- data.frame()

	printlog ('Opening file', fn )
	ncid <- open.ncdf( fn )
	print( ncid )

	printlog( "Getting lat/lon data" )  
	lat <- get.var.ncdf( ncid, LAT_NAME )
	lon <- get.var.ncdf( ncid, LON_NAME )
	# lat and lon can be returned in two different ways
	if( is.matrix( lat ) ) {
		latlon <- data.frame( Var1=melt( lon )$value, Var2=melt( lat )$value )
	} else {
		latlon <- expand.grid( lon, lat )    
	}  
	
	tf <- tempfile()
	printlog( "Using tempfile", tf )

	dnames <- read_dimension_names( ncid, variable )
	printlog( "Dimensions names for", variable, "=", dnames )

	levels_to_avg <- determine_levels( ncid, dnames, LEVEL_DEPTH )

	latindex <- which( dnames==LAT_NAME | dnames==LAT_NAME2 | dnames==LAT_NAME3 | dnames==LAT_NAME4 )
	lonindex <- which( dnames==LON_NAME | dnames==LON_NAME2| dnames==LON_NAME3 | dnames==LON_NAME4)
	levindex <- which( dnames==LEVEL_NAME )
	timeindex <- which( dnames==TIME_NAME )
	printlog( "latindex lonindex levindex timeindex" )
	printlog( latindex, lonindex, levindex, timeindex )
	startdata <- rep( 1, length( dnames ) )
	countdata <- rep( -1, length( dnames ) )
	countdata[ timeindex ] <- 1 # always reading only 1 time slice at a time
	countdata[ levindex ] <- 1  # always reading only 1 level at a time

	first_time <- T
    year_range <- beginyear:endyear
    for( yearindex in 1:length( year_range ) ) {
        year <- year_range[ yearindex ]
        inrange <- year %in% YEAR_RANGE1 | year %in% YEAR_RANGE2
        if( !inrange ) {
           # printlog( "-- skipping year", year )
            next
        }
        printlog( "-- processing year index", yearindex, year )
        
        for( month in 1:12 ) {		# assume we're not processing any fractional years
            #printlog( "--   reading month", month )
            startdata[ timeindex ] <- ( year_range[ yearindex ]-beginyear ) * MONTHS_PER_YEAR + month
            d <- read_timepoint( ncid, variable, levels_to_avg, startdata, countdata, levindex )
            d_m <- melt( d )
            d_m$year <- year_range[ yearindex ]
            d_m$lat <- latlon$Var2
            d_m$lon <- latlon$Var1
            
            d_m$month <- month
            
           write.table( d_m, file=tf, row.names=F, col.names=first_time, sep=",", append = !first_time )
            first_time <- F
        } # for month
    } # for year

	if( file.exists( tf ) ) {
	    printlog( "Reading tempfile back into results..." )
	    results <- read.csv( tf )
	    printdims( results )
	    printlog( "Size =", format( object.size( results ), units = "Mb" ) )
	    # printlog( "Removing NA values and rounding..." )
	    results <- subset( results, !is.na( value ) )
	    # results$value <- round( results$value, 6 )
	    results$units <- att.get.ncdf (ncid, variable,  "units" )$value
	} else {
	    results <- NULL
	}
    
	close.ncdf( ncid )
    
    return( results ) 
} # process_data

# -----------------------------------------------------------------------------
# Parse a CMIP5 filename
parse_filename <- function( fn ) {
} # parse_filename

# -----------------------------------------------------------------------------
# Process a single file, parsing information from its name and calling process_data
# allow_historical is a flag indiating whether we should process this file if it's historical
# normally no: we skip historical files
# but when recursing, can't find years in a scenario, allow and return data to caller (us)
process_file <- function( fn, tf ) {
	filedata <- strsplit( basename( fn ), "_" )[[ 1 ]]
	printlog( "------------------------" )
	printlog( "File:", fn )
	variable <- filedata[ 1 ]
	model <- filedata[ 3 ]
	scenario <- filedata[ 4 ]
	ensemble <- filedata[ 5 ]

    if( tolower( scenario ) != tolower( SCENARIO ) ) {
        if( tolower( scenario ) != "historical" ) {
            printlog( "File is not", SCENARIO, "or historical--skipping" )
            return( NULL )
        }
    }
	printlog( variable, model, scenario, ensemble )
	filename <- strsplit( filedata[ 6 ], ".", fixed=T )[[ 1 ]][ 1 ]
	timedata <- strsplit( filename, "-" )[[ 1 ]]
	begintime <- timedata[ 1 ]
	endtime <- timedata[ 2 ]

	if( nchar( begintime ) != nchar( endtime ) | ( nchar( begintime ) != 4 & nchar( begintime ) != 6 ) ) {
		printlog( "*** Uh oh! Something's wrong--date not 4 or 6 digits. Skipping ***" )
		invfile( fn, status="Date not 4 or 6 digits" )
		warning( paste( "Couldn't parse filename", i ) )
		return( NULL )
	}

	beginyear <- as.numeric( substr( begintime, 1, 4 ) )
	endyear <- as.numeric( substr( endtime, 1, 4 ) )
	beginmonth <- 0
	endmonth <- 0

	if( nchar( begintime )==6 ) {
		datafreq <- DATAFREQ_MONTHLY
		beginmonth <- as.numeric( substr( begintime, 5, 6 ) )
		endmonth <- as.numeric( substr( endtime, 5, 6 ) )
	} else {
		datafreq <- DATAFREQ_ANNUAL
	}

	printlog( "This appears to be", datafreq, "data" )
	printlog( beginyear, beginmonth, endyear, endmonth )

	results <- process_data( fn, variable, beginyear, endyear, datafreq )

    
    if( !is.null( results ) ) {
        results$model <- model
        #    results$scenario <- scenario
        results$ensemble <- ensemble
        printlog( "Writing data to tempfile..." )
        first <- !file.exists( tf )
        write.table( results, file=tf, sep=",", row.names=F, col.names=first, append=!first )        
    }
} # process_file

# -----------------------------------------------------------------------------
# Process a particular variable, in a (recursive) directory of files
process_directory <- function( model, variable, dir=INPUT_DIR, pattern="*.nc$" ) {

	invfile( newfile=T )	# erase the skip log and start a new one

	printlog( SEPARATOR )
	printlog( "Welcome to process_directory: doing", model, variable )
    
    # filter file list to just the variable and model we're processing
	files <- list.files( dir, pattern=pattern, recursive=TRUE )
	files <- files[ grepl( paste0( "^", variable ), basename(files) ) ]
	files <- files[ grepl( paste0( model ), files ) ]
	
	printlog( length( files ), "files to process" )
    
    if( any( duplicated( basename( files ) ) ) ) {
     printlog( "WARNING...Duplicate files found:", files[ duplicated( basename( files ) ) ] )
     stop( 'Fix this!' )
    }
    
    tf <- tempfile()
    printlog( "Tempfile:", tf )

	for( i in files ) {
		process_file( paste( dir, i, sep="/" ), tf )
	}
    
	printlog( SEPARATOR )
	if( !file.exists( tf ) ) {
        printlog( "No data written to tempfile! Exiting" )
        return()
    }
    
	printlog( "Reading tempfile for", model, variable, "back in..." )
    results <- read.csv( tf )
    printdims( results )
    results$range <- NA
    results[ results$year %in% YEAR_RANGE1, "range" ] <- 1
	results[ results$year %in% YEAR_RANGE2, "range" ] <- 2


    printlog( "Checking data..." )
    results1 <- subset( results, year %in% YEAR_RANGE1 )
    results2 <- subset( results, year %in% YEAR_RANGE2 )
	printdims( results1 )
	printdims( results2 )
    
	if( nrow( results1 ) != nrow( results2 ) ) {
	    printlog( "*** Uh oh! Something's wrong--Unequal data sets for year ranges ***" )
	    invfile( variable, status="Unequal data sets for year ranges" )
        printlog("Data are in", tf)
        return()
	}
	

    # Calculate a single monthly mean across all years for each grid point
	printlog( "Aggregating years...." )
	print( system.time( 
{ 	# base R: aggregate
    #results_agg <- aggregate( value~month+lat+lon, data=results, mean )	 # aggregate faster than ddply 
    
    # plyr
    #results_agg <- ddply( results, .( lat, lon, month ), summarise, value=mean( value ), nyears=length( lat ), .progress="text" )
    
    # dplyr
    results_grouped <- group_by( results, model, ensemble, range, month, lat, lon, units )
    results_agg <- dplyr::summarise( results_grouped, value=mean( value ), year=mean( year ) )
}
	) )

    printdims( results_agg )

    # compute differences between the two year ranges
    results1 <- subset( results_agg, range==1 )
    results2 <- subset( results_agg, range==2 )
    results3 <- results1
    results3$value <- results2$value - results1$value
    results3$year <- length( YEAR_RANGE2 )
#    results3$range <- 3
    final_results <- rbind( results_agg, results3 )
    final_results$range <- NULL

    for( m in unique( final_results$model ) ) {
        outfnm <- paste0( OUTPUT_DIR, variable, "_", SCENARIO, "_", m, ".csv" )
        d <- subset( final_results, model==m )
        d$model <- NULL
        printlog( "Writing", outfnm )
        write.csv( d, file=outfnm, row.names=F )
    }
#    invfile( fn, model, scenario, ensemble, units )

} # process_directory

# =============================================================================
# Main

if( !file.exists( OUTPUT_DIR ) ) {
	printlog( "Creating", OUTPUT_DIR )
	dir.create( OUTPUT_DIR )
}
if( !file.exists( LOG_DIR ) ) {
	printlog( "Creating", LOG_DIR )
	dir.create( LOG_DIR )
}

sink( paste0( LOG_DIR, SCRIPTNAME, ".txt" ), split=T )

printlog( "Welcome to", SCRIPTNAME )

printlog( "We are processing years", range( YEAR_RANGE1 ), "and", range( YEAR_RANGE2 ) )

loadlibs( c( "ncdf", "reshape2", "ggplot2", "plyr", "dplyr" ) )
theme_set( theme_bw() )

for( m in MODELS ) {
    for( v in VARIABLES ) {
        process_directory( m, v )
    }
}

printlog( "All done." )
sink()
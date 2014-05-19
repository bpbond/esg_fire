# Script to process ESG files into CSV annual summary format
# B.B.-L. and C.H. for Yannick fire project
# March 2014
# combines both land and ocean data into one script with latitude and longitude bands

# Important variable definitions, esp. data source & destination
SCRIPTNAME		<- "esg_data_fire.R"
INPUT_DIR		<- "sampledata/Corinne"
OUTPUT_DIR		<- "outputs/"
HIST_DIR        <- "../historical/"  # relative to location of file being processed
LOG_DIR			<- "logs/"
SEPARATOR		<- "-------------------"

# To use:
#	1. Call process_directory or process_file as appropriate

YEAR_RANGE1			<- 1996:1999
YEAR_RANGE2			<- 2066:2069	# these should be same length

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
}

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
}

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
}

# -----------------------------------------------------------------------------
# The workhorse: read a file's data, fill data frame, write out
process_data <- function( fn, variable, beginyear, beginmonth, endyear, endmonth, datafreq, year_range ) {

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

	printlog( "Year range looking for is", year_range )
	
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

	for( yearindex in 1:length( year_range ) ) {
		printlog( "-- processing year index", yearindex )

		for( month in 1:12 ) {		# assume we're not processing any fractional years
			printlog( "--   reading month", month )
			startdata[ timeindex ] <- ( MONTHS_PER_YEAR-beginmonth+1 ) + ( year_range[ yearindex ]-beginyear-1 ) * MONTHS_PER_YEAR + month

			d <- read_timepoint( ncid, variable, levels_to_avg, startdata, countdata, levindex )
			d_m <- melt( d )
			d_m$year <- year_range[ yearindex ]
			d_m$lat <- latlon$Var2
			d_m$lon <- latlon$Var1

			d_m$month <- month

			write.table( d_m, file=tf, row.names=F, col.names=( yearindex==1 ), sep=",", append =! ( yearindex==1 ) )
		} # for month
	} # for year

	printlog( "Reading tempfile back into results..." )
	results <- read.csv( tf )
	printdims( results )
	printlog( "Size =", format( object.size( results ), units = "Mb" ) )
	printlog( "Removing NA values and rounding..." )
	results <- subset( results, !is.na( value ) )
	results$value <- round( results$value, 6 )
	results$X1 <- NULL
	results$X2 <- NULL
	results$units <- att.get.ncdf (ncid, variable,  "units" )$value
	results$variable <- variable

	close.ncdf( ncid )

	return( results ) 
}

# -----------------------------------------------------------------------------
# Parse a CMIP5 filename
parse_filename <- function( fn ) {
}

# -----------------------------------------------------------------------------
# Process a single file, parsing information from its name and calling process_data
# allow_historical is a flag indiating whether we should process this file if it's historical
# normally no: we skip historical files
# but when recursing, can't find years in a scenario, allow and return data to caller (us)
process_file <- function( fn, skip_existing=FALSE, allow_historical=F ) {
	filedata <- strsplit( basename( fn ), "_" )[[ 1 ]]
	printlog( "------------------------" )
	printlog( "File:", fn )
	variable <- filedata[ 1 ]
	model <- filedata[ 3 ]
	scenario <- filedata[ 4 ]
	ensemble <- filedata[ 5 ]
	outfn <- paste0( OUTPUT_DIR, basename( fn ), ".csv" )
	
	if( scenario=="historical" & !allow_historical ) {
		printlog( "Skipping this historical file" )
		invfile( fn, status="Historical file" )
		return( NULL )
	}
	
	if (file.exists( outfn ) & skip_existing ) {
		printlog ("Skipping file", fn )
		invfile( fn, "Output filename already exists" )
		return( NULL )
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

	results1 <- NULL  
	results2 <- NULL  
	if( all( YEAR_RANGE1 %in% beginyear:endyear ) ) {
		results1 <- process_data( fn, variable, beginyear, beginmonth, endyear, endmonth, datafreq, YEAR_RANGE1 )
	} else {
		printlog( "Skipping YEAR_RANGE1 - dates not included in this file" )
		
		# ...but these years might be in a 'historical' file
		# if so, recurse to get those data
		otherfn <- paste( variable, filedata[ 2 ], model, "historical", ensemble, "*.*", sep="_" )
        otherdir <- paste0( dirname( fn ), "/", HIST_DIR)
        printlog ("Looking in folder", otherdir)
        
		filelist <- list.files( otherdir, otherfn )
		if( length( filelist )==1 & !allow_historical ) {
			otherfn <- filelist[ 1 ]
			printlog( "Possible alternate file exists:", otherfn )
#			readline()
			printlog( "Shifting to historical file" )
			results1 <- process_file( paste0( otherdir, "/", otherfn ), skip_existing, allow_historical=T )
			printlog( "We're back" )
			print( summary( results1 ) )
		} else {
			printlog ("no alternative file found")
		}
	}
	
	if( all( YEAR_RANGE2 %in% beginyear:endyear ) ) {
		results2 <- process_data( fn, variable, beginyear, beginmonth, endyear, endmonth, datafreq, YEAR_RANGE2 )
	} else {
		printlog( "Skipping YEAR_RANGE2 - dates not included in this file" )
	}
	
	results <- rbind( results1, results2 )

	if( allow_historical ) {		# bail - we don't want to write anything out
		printlog( "Returning data to caller without writing..." )
		return( results )
	}
	
	# compute differences between results1 and results2, if both available
	if( !is.null( results1 ) & !is.null( results2 ) ) {
		results3 <- results1
		stopifnot( nrow( results1 )==nrow( results2 ) )
		# we're assuming that the two data frames are structured identically
		results3$value <- results2$value - results1$value
		results3$year <- results3$year - min( results3$year )
		results <- rbind( results, results3 )
	}
	
	if( !is.null( results ) ) {
		units <- results[ 1, "units" ]
		results$units <- NULL		# these will go to log file instead
		results$variable <- NULL		# these will go to log file instead

		printlog( "Writing output file", outfn, "..." )
		write.csv( results, file=outfn, row.names=F )
		invfile( fn, model, scenario, ensemble, units )
	} else {
		printlog( "NULL results! Not writing any output" )
		invfile( fn, status="NULL results back from process_data" )
	} # if
	return( NULL )
}

# -----------------------------------------------------------------------------
# Process a whole directory of files
process_directory <- function( dir, pattern="*.nc$" ) {

	invfile( newfile=T )	# erase the skip log and start a new one

	printlog( SEPARATOR )
	printlog( "Welcome to process_directory_step1" )
	files <- list.files( dir, pattern=pattern )
	printlog( length( files ), "files to process" )

	for( i in files ) {
		process_file( paste( dir, i, sep="/" ) )
	}
}

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

loadlibs( c( "ncdf", "reshape", "ggplot2", "plyr" ) )
theme_set( theme_bw() )

process_directory( INPUT_DIR )

printlog( "All done." )
sink()
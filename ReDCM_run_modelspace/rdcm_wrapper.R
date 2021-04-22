#!/usr/bin/env Rscript


## -----------------------------------------------------------------------------------------
## set environmental variables
.libPaths(rev(dir('~/R/x86_64-pc-linux-gnu-library', full.names=TRUE)))
cat( 'R package path: ' )
cat( .libPaths(), '\n', sep=':')


suppressPackageStartupMessages( library("methods") )
suppressPackageStartupMessages( library("optparse") )
suppressPackageStartupMessages( library("ReDCM") )


## -----------------------------------------------------------------------------------------
## functions (source absolute paths)

thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

w.script.dir = dirname( thisFile() )
source( str_c(w.script.dir, '/rdcm_wrapper_util.R') )
cat('script.dir (wrapper): ', w.script.dir, '\n')


## -----------------------------------------------------------------------------------------
## globals

#data.path = "/data"
#accepted.sereieses.file = "/data/etc/accpeted-sereieses.csv"




## =====================================================================================
## options
option_list = list(
    
  make_option(c("-a", "--amatrix"), type="character", default=NULL,
              help="A matrix file", metavar="character"),
  
  make_option(c("-b", "--bmatrix"), type="character", default=NULL,
              help="B matrix file", metavar="character"),
  
  make_option(c("-c", "--cmatrix"), type="character", default=NULL,
              help="C matrix file", metavar="character"),
  
  make_option(c("-d", "--dmatrix"), type="character", default=NULL,
              help="D matrix file", metavar="character"),
  
  make_option(c("-x", "--xsubspace"), type="integer", default=1,
              help="model space interval begin", metavar="number"),
  
  make_option(c("-y", "--ysubspace"), type="integer", default=NULL,
              help="model space interval end", metavar="number"),
  
  make_option(c("-n", "--noffset"), type="integer", default=0,
              help="offset of output DCM name", metavar="number"),
  
  make_option(c("-o", "--output_dir"), type="character", default='./out',
              help="output directory [default = %default]", metavar="character")
  
);

## set parser
parser    = OptionParser(usage = "%prog DCM_template [options]", option_list=option_list);

## check arguments
args = commandArgs(trailingOnly = TRUE)

if ( length(args) == 0 )
{
  print_help( parser )
  stop("this script requires DCM.mat template as argument")
}

arguments = parse_args(parser, args = args, positional_arguments = 1)

w.options = arguments$options
w.dcmarg = arguments$args

#cat(options[[1]], '\n')
#cat(path[[1]], '\n')

## =====================================================================================
## MAIN


# 1. check input filenames
# -------------------------------------

w.sessID = NULL

if ( is.abs.path(w.options$o) )
{
  dir.create( w.options$o, recursive=TRUE, showWarnings=FALSE )
  cat('outpath: ', w.options$o, '\n')
} else
{
  w.options$o = str_c( w.script.dir, w.options$o, sep='/')
  dir.create( w.options$o, recursive=TRUE, showWarnings=FALSE )
  cat('outpath: ', w.options$o, '\n')
}

if ( !is.abs.path(w.dcmarg[[1]]) )
{
  w.dcmtemp = str_c( w.script.dir, w.dcmarg[[1]], sep='/')
  cat('dcm template: ', w.dcmtemp, '\n')
} else
{
  w.dcmtemp = w.dcmarg[[1]]
  cat('dcm template: ', w.dcmtemp, '\n')
}



# 2. read DCM structure
# -------------------------------------

w.DCM = readMat(w.dcmtemp)
w.DCM = w.DCM$DCM[,,1]
w.DCM.options = w.DCM$options[,,1]
w.DCM.n = length(w.DCM$Y[,,1]$name)
w.DCM.uN = length(w.DCM$U[,,1]$name)


# 3. read options - check matrix filenames and read matrices
# -------------------------------------

w.matrices = list(NULL, NULL, NULL, NULL)
w.havemx = NULL
w.mxlen = c(0,0,0,0)

if ( !is.null(w.options$a) )
{
  if ( !is.abs.path(w.options$a) )
  {
    w.options$a = str_c( w.script.dir, w.options$a, sep='/')
    cat('A matrix: ', w.options$a, '\n')
  }
  w.matrices[[1]] = read.table( w.options$a, colClasses='character')
  w.havemx = c(w.havemx, 1)
  w.mxlen[1] = length(w.matrices[[1]][,1])
}
if ( !is.null(w.options$b) )
{
  if ( !is.abs.path(w.options$b) )
  {
    w.options$b = str_c( w.script.dir, w.options$b, sep='/')
    cat('B matrix: ', w.options$b, '\n')
  }
  w.matrices[[2]] = read.table( w.options$b, colClasses='character')
  w.havemx = c(w.havemx, 2)
  w.mxlen[2] = length(w.matrices[[2]][,1])
}
if ( !is.null(w.options$c) )
{
  if ( !is.abs.path(w.options$c) )
  {
    w.options$c = str_c( w.script.dir, w.options$c, sep='/')
    cat('C matrix: ', w.options$c, '\n')
  }
  w.matrices[[3]] = read.table( w.options$c, colClasses='character')
  w.havemx = c(w.havemx, 3)
  w.mxlen[3] = length(w.matrices[[3]][,1])
}
if ( !is.null(w.options$d) )
{
  if ( !is.abs.path(w.options$d) )
  {
    w.options$d = str_c( w.script.dir, w.options$d, sep='/')
    cat('D matrix: ', w.options$d, '\n')
  }
  w.matrices[[4]] = read.table( w.options$d, colClasses='character')
  w.havemx = c(w.havemx, 4)
  w.mxlen[4] = length(w.matrices[[4]][,1])
}

#cat(w.mxlen, '\n')
#cat(w.havemx, '\n')
#cat('matrix sizes: ', unique( w.mxlen[w.havemx] ), '\n')
#stop("emergency abort")


# 4. prepare and estimate model space
# -------------------------------------

if ( !is.null(w.havemx) )
{
  # estimate model subspace
  #if ( length(unique( w.mxlen[w.havemx] )) > 1 )
  #{
  #  stop("matrix files contains different number of rows")
  #}
  
  # check model space interval
  if ( is.null(w.options$y) )
  {
    w.options$y = w.options$x
  }
  
  if ( w.options$y > unique( w.mxlen[w.havemx] ) )
  {
    w.options$y = unique( w.mxlen[w.havemx] )
  }
  if ( w.options$x > unique( w.mxlen[w.havemx] ) )
  {
    w.options$x = unique( w.mxlen[w.havemx] )
  }
    
  #if ( (w.options$y > unique( w.mxlen[w.havemx] )) || (w.options$x < 0) || (w.options$y < 0) || (w.options$y < w.options$x) )
  if ( (w.options$x < 0) || (w.options$y < 0) || (w.options$y < w.options$x) )
  {
    stop("wrong model space interval")
  }
  
  cat('Model subspace:: [', w.options$x, ',', w.options$y, ']\n\n')
  cat('preA: ', c(w.DCM$a), '\n', sep='')
  cat('preB: ', c(w.DCM$b), '\n', sep='')
  cat('preC: ', c(w.DCM$c), '\n', sep='')
  cat('preD: ', c(w.DCM$d), sep='')
  #stop("emergency abort")
  
  # cycle through model subspace
  w.modelList = seq(w.options$x, w.options$y)
  for ( w.m in w.modelList )
  {

    w.fn = str_c(w.options$o, '/DCM_', w.m+w.options$n)
    if (file.exists(w.fn))
    {
      cat(str_c('\n', w.options$o, '/DCM_', w.m+w.options$n, ' already exists'))
      next
    }

    cat('\n\nCurrent model: ', w.m, '\n')

    for ( w.i in w.havemx )
    {
      # swap A matrix
      if ( w.i == 1 )
      {
        w.tmp = as.numeric(str_split(w.matrices[[w.i]][w.m,], pattern='')[[1]])
        w.tmplen = sqrt(length(w.tmp))
        if ( w.tmplen != w.DCM.n )
        {
          stop("wrong A matrix size")
        }
        else
        {
          w.DCM$a = array( w.tmp, dim=c(w.DCM.n, w.DCM.n) )
        }
      }
      # swap B matrix
      if ( w.i == 2 )
      {
        w.tmp = as.numeric(str_split(w.matrices[[w.i]][w.m,], pattern='')[[1]])
        w.tmplen = sqrt(length(w.tmp) / w.DCM.uN)
        if ( w.tmplen != w.DCM.n )
        {
          stop("wrong B matrix size")
        }
        else
        {
          w.DCM$b = array( w.tmp, dim=c(w.DCM.n, w.DCM.n, w.DCM.uN) )
        }
      }
      # swap C matrix
      if ( w.i == 3 )
      {
        w.tmp = as.numeric(str_split(w.matrices[[w.i]][w.m,], pattern='')[[1]])
        w.tmplen = length(w.tmp) / w.DCM.uN
        if ( w.tmplen != w.DCM.n )
        {
          stop("wrong C matrix size")
        }
        else
        {
          w.DCM$c = array( w.tmp, dim=c(w.DCM.n, w.DCM.uN) )
        }
      }
      # swap D matrix
      if ( w.i == 4 )
      {
        w.tmp = as.numeric(str_split(w.matrices[[w.i]][w.m,], pattern='')[[1]])
        w.tmplen = sqrt(length(w.tmp) / w.DCM.n)
        if ( w.tmplen != w.DCM.n )
        {
          stop("wrong D matrix size")
        }
        else
        {
          w.DCM$d = array( w.tmp, dim=c(w.DCM.n, w.DCM.n, w.DCM.n) )
        }
      }
    }
    
    cat('   A: ', c(w.DCM$a), '\n', sep='')
    cat('   B: ', c(w.DCM$b), '\n', sep='')
    cat('   C: ', c(w.DCM$c), '\n', sep='')
    cat('   D: ', c(w.DCM$d), '\n\n', sep='')
    #stop("emergency abort")
    
    #estimate current DCM model
    
    w.DCMe = ReDCM_prepare_dcm( w.DCM )
    
    w.ptm = proc.time()
    w.DCMe = ReDCM_estimate( w.DCMe )
    w.ptm = proc.time() - w.ptm
    w.DCMe = setupEstimates( w.DCMe, tt=w.ptm[3] )
    
    w.fn = str_c(w.options$o, '/DCM_', w.m+w.options$n)
    DCMe = w.DCMe
    save(DCMe, file=w.fn, ascii=FALSE)
    rm( DCMe )
    rm( w.DCMe )
    
    #stop("emergency abort")
    #rdcm_print_output_csv( w.DCMe, w.options$o )
    
  }
  
} else
{
  # estimate DCM.mat template
  w.ptm = proc.time()
  DCMe = rdcm_estimate( w.dcmtemp )
  w.ptm = proc.time() - w.ptm
  
  w.DCMe = setupEstimates( w.DCMe, tt=w.ptm[3] )
  
  w.fn = str_c(w.options$o, '/DCM_', w.m+w.options$n)
  DCMe = w.DCMe
  save(DCMe, file=w.fn, ascii=FALSE)
  rm( DCMe )
  rm( w.DCMe )
  
  #rdcm_print_output_csv( w.DCMe, w.options$o )
}

cat('\n')

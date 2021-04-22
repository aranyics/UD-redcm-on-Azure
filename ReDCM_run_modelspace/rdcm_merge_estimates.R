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

setwd(w.script.dir)

cat('\n')

## -----------------------------------------------------------------------------------------
## globals




## =====================================================================================
## options
option_list = list(
  
  make_option(c("-x", "--xsubspace"), type="integer", default=1,
              help="model space interval begin", metavar="number"),
  
  make_option(c("-y", "--ysubspace"), type="integer", default=NULL,
              help="model space interval end", metavar="number"),
  
  make_option(c("-p", "--prefix"), type="character", default='DCM',
              help="prefix of estimated DCM stucture filenames [default = %default]", metavar="character"),
  
  make_option(c("-o", "--output"), type="character", default='csv/',
              help="output path of merged .csv files [default = %default]", metavar="character")
  
);

## set parser
parser    = OptionParser(usage = "%prog DCM_dir [options]", option_list=option_list);

## check arguments
args = commandArgs(trailingOnly = TRUE)

if ( length(args) == 0 )
{
  print_help( parser )
  stop("this script requires DCM output directory as argument")
}

arguments = parse_args(parser, args = args, positional_arguments = 1)

w.options = arguments$options
w.dcmarg = arguments$args

#cat(options[[1]], '\n')
#cat(path[[1]], '\n')

## =====================================================================================
## MAIN

# 1. check input paths and filenames
# -------------------------------------

w.sessID = NULL
w.dcm.dir = NULL

if ( !is.abs.path(w.dcmarg[[1]]) )
{
  w.dcm.dir = str_c( w.script.dir, w.dcmarg[[1]], sep='/')
  cat('dcm output dir: ', w.dcm.dir, '\n')
} else
{
  w.dcm.dir = w.dcmarg[[1]]
  cat('dcm output dir: ', w.dcm.dir, '\n')
}

w.dcm.files = dir( w.dcm.dir )[which(grepl( dir( w.dcm.dir ), pattern=str_c('^',w.options$p)))]
if ( length(w.dcm.files) < 1 )
{
  stop("no DCM structure files found")
}

# 2. check model space interval
# -------------------------------------

if ( (w.options$x < 0) )
{
  w.options$x = 0
}

if ( is.null(w.options$y) )
{
  w.dcm.exclude = str_c(w.options$p, seq(0, w.options$x-1), sep='_')
  w.tmp = match( w.dcm.exclude, w.dcm.files )
  if ( any(!is.na(w.tmp)) )
  {
    w.dcm.files = w.dcm.files[ -w.tmp[which( !is.na(w.tmp) )] ]
  }
} else
{
  if ( (w.options$y < w.options$x) )
  {
    w.options$y = w.options$x
  }
  
  w.dcm.include = str_c(w.options$p, seq(w.options$x, w.options$y), sep='_')
  w.tmp = match( w.dcm.include, w.dcm.files )
  w.dcm.files = w.dcm.files[ w.tmp[which( !is.na(w.tmp) )] ]  
}

if ( is.null( w.options$y ) )
{
  cat('Model subspace:: [', w.options$x, ', end ]\n\n')
} else
{
  cat('Model subspace:: [', w.options$x, ',', w.options$y, ']\n\n')
}
  
#stop("emergency abort")



# 3. create merged csv from estimated DCMs
# -------------------------------------

w.grp = 'na'
w.id = 'na'
w.condition = 'r01'
w.task = 'vdmodel'
w.rep = '1'
w.dcm = 'd1l'
w.F = 0
w.tt = 0
w.iter = 0
w.matrix = 'mx'

setwd(w.dcm.dir)
cat( 'Creating merged results -> ', w.dcm.dir, '/', w.options$o, '/', w.options$p, 'e_*_', w.options$x, '-', w.options$y, '.csv from:\n', sep='' )
dir.create(w.options$o, showWarnings=FALSE)
w.fl.head=FALSE

# print matrices by models
for ( w.m in w.dcm.files )
{
  
  cat( w.m, '\n' )
  load( str_c(w.dcm.dir, w.m, sep='/') )
  
  if ( exists('w.DCMe') )
  {
    DCM = w.DCMe
    rm(w.DCMe)
  } else if ( exists('DCMe') )
  {
    DCM = DCMe
    rm(DCMe)
  } else
  {
    stop("couldn't load DCM structure")
  }
  
  # load data from DCM struct
  w.F = DCM@Fe
  w.tt = DCM@tt
  w.iter = DCM@k
  
  # print headers
  if ( !w.fl.head )
  {
    w.head = t( c('grp', 'id', 'condition', 'task', 'repetition', 'dcm', 'F', 'time', 'iter', 'matrix', colnames.cross(DCM@Y@name)) )
    w.fn = str_c( w.options$o, '/DCMe_A_', w.options$x, '-', w.options$y, '.csv' )
    write.table( w.head, file=w.fn, sep=',', row.names=FALSE, col.names=FALSE )
    w.fn = str_c( w.options$o, '/DCMe_B_', w.options$x, '-', w.options$y, '.csv' )
    if ( length(DCM@Ep@D) > 0 )
    {
      write.table( w.head, file=w.fn, sep=',', row.names=FALSE, col.names=FALSE )
      w.fn = str_c( w.options$o, '/DCMe_D_', w.options$x, '-', w.options$y, '.csv' )
    }
    write.table( w.head, file=w.fn, sep=',', row.names=FALSE, col.names=FALSE )
    w.head = t( c('grp', 'id', 'condition', 'task', 'repetition', 'dcm', 'F', 'time', 'iter', 'matrix', DCM@Y@name) )
    w.fn = str_c( w.options$o, '/DCMe_C_', w.options$x, '-', w.options$y, '.csv' )
    write.table( w.head, file=w.fn, sep=',', row.names=FALSE, col.names=FALSE )
    w.head = t( c('grp', 'id', 'condition', 'task', 'repetition', 'dcm', 'F', 'time', 'iter', 'matrix', 
                  c( paste('D', DCM@Y@name, sep='.'), paste('T', DCM@Y@name, sep='.'), 'E' )) )
    w.fn = str_c( w.options$o, '/DCMe_H_', w.options$x, '-', w.options$y, '.csv' )
    write.table( w.head, file=w.fn, sep=',', row.names=FALSE, col.names=FALSE )
    w.fl.head=TRUE
  }
  
  # print A and pA matrices
  w.fn = str_c( w.options$o, '/DCMe_A_', w.options$x, '-', w.options$y, '.csv' )
  write.table( data.frame( grp=w.grp,
                           id=str_split_fixed(w.m, pattern = "_", 2)[2],
                           condition=w.condition,
                           task=w.task,
                           repetition=w.rep,
                           dcm=w.dcm,
                           F=w.F,
                           time=w.tt,
                           iter=w.iter,
                           matrix='A',
                           t( c(DCM@Ep@A) ) ),
               file=w.fn,
               sep=',',
               col.names=FALSE,
               row.names=FALSE,
               append=TRUE )
  write.table( data.frame( grp=w.grp,
                           id=str_split_fixed(w.m, pattern = "_", 2)[2],
                           condition=w.condition,
                           task=w.task,
                           repetition=w.rep,
                           dcm=w.dcm,
                           F=w.F,
                           time=w.tt,
                           iter=w.iter,
                           matrix='pA',
                           t( c(DCM@Pp@A) ) ),
               file=w.fn,
               sep=',',
               col.names=FALSE,
               row.names=FALSE,
               append=TRUE )
  
  # print Bj and pBj matrices
  for (w.j in 1 : DCM@M[[1]]@m)
  {
    w.fn = str_c( w.options$o, '/DCMe_B_', w.options$x, '-', w.options$y, '.csv' )
    write.table( data.frame( grp=w.grp,
                             id=str_split_fixed(w.m, pattern = "_", 2)[2],
                             condition=w.condition,
                             task=w.task,
                             repetition=w.rep,
                             dcm=w.dcm,
                             F=w.F,
                             time=w.tt,
                             iter=w.iter,
                             matrix=str_c('B', w.j, sep=''),
                             t( c(DCM@Ep@B[,,w.j]) ) ),
                 file=w.fn,
                 sep=',',
                 col.names=FALSE,
                 row.names=FALSE,
                 append=TRUE )
    write.table( data.frame( grp=w.grp,
                             id=str_split_fixed(w.m, pattern = "_", 2)[2],
                             condition=w.condition,
                             task=w.task,
                             repetition=w.rep,
                             dcm=w.dcm,
                             F=w.F,
                             time=w.tt,
                             iter=w.iter,
                             matrix=str_c('pB', w.j, sep=''),
                             t( c(DCM@Pp@B[,,w.j]) ) ),
                 file=w.fn,
                 sep=',',
                 col.names=FALSE,
                 row.names=FALSE,
                 append=TRUE )  
  
    # print C and pC matrices
    w.fn = str_c( w.options$o, '/DCMe_C_', w.options$x, '-', w.options$y, '.csv' )
    write.table( data.frame( grp=w.grp,
                             id=str_split_fixed(w.m, pattern = "_", 2)[2],
                             condition=w.condition,
                             task=w.task,
                             repetition=w.rep,
                             dcm=w.dcm,
                             F=w.F,
                             time=w.tt,
                             iter=w.iter,
                             matrix=str_c('C', w.j, sep=''),
                             t( c(DCM@Ep@C[,w.j]) ) ),
                 file=w.fn,
                 sep=',',
                 col.names=FALSE,
                 row.names=FALSE,
                 append=TRUE )
    write.table( data.frame( grp=w.grp,
                             id=str_split_fixed(w.m, pattern = "_", 2)[2],
                             condition=w.condition,
                             task=w.task,
                             repetition=w.rep,
                             dcm=w.dcm,
                             F=w.F,
                             time=w.tt,
                             iter=w.iter,
                             matrix=str_c('pC', w.j, sep=''),
                             t( c(DCM@Pp@C[,w.j]) ) ),
                 file=w.fn,
                 sep=',',
                 col.names=FALSE,
                 row.names=FALSE,
                 append=TRUE )
  }
  
  # print Di and pDi matrices
  if ( length(DCM@Ep@D) > 0 )
  {
    for (w.i in 1 : DCM@M[[1]]@l)
    {
      w.fn = str_c( w.options$o, '/DCMe_D_', w.options$x, '-', w.options$y, '.csv' )
      write.table( data.frame( grp=w.grp,
                               id=str_split_fixed(w.m, pattern = "_", 2)[2],
                               condition=w.condition,
                               task=w.task,
                               repetition=w.rep,
                               dcm=w.dcm,
                               F=w.F,
                               time=w.tt,
                               iter=w.iter,
                               matrix=str_c('D', w.i, sep=''),
                               t( c(DCM@Ep@B[,,w.i]) ) ),
                   file=w.fn,
                   sep=',',
                   col.names=FALSE,
                   row.names=FALSE,
                   append=TRUE )
      write.table( data.frame( grp=w.grp,
                               id=str_split_fixed(w.m, pattern = "_", 2)[2],
                               condition=w.condition,
                               task=w.task,
                               repetition=w.rep,
                               dcm=w.dcm,
                               F=w.F,
                               time=w.tt,
                               iter=w.iter,
                               matrix=str_c('pD', w.i, sep=''),
                               t( c(DCM@Pp@B[,,w.i]) ) ),
                   file=w.fn,
                   sep=',',
                   col.names=FALSE,
                   row.names=FALSE,
                   append=TRUE ) 
    }
  }
  
  # print H and pH hemodynamic parameters
  w.fn = str_c( w.options$o, '/DCMe_H_', w.options$x, '-', w.options$y, '.csv' )
  write.table( data.frame( grp=w.grp,
                           id=str_split_fixed(w.m, pattern = "_", 2)[2],
                           condition=w.condition,
                           task=w.task,
                           repetition=w.rep,
                           dcm=w.dcm,
                           F=w.F,
                           time=w.tt,
                           iter=w.iter,
                           matrix='H',
                           t( c( t(c(DCM@Ep@decay)),t(c(DCM@Ep@transit)),DCM@Ep@epsilon ) )  ),
               file=w.fn,
               sep=',',
               col.names=FALSE,
               row.names=FALSE,
               append=TRUE )
  write.table( data.frame( grp=w.grp,
                           id=str_split_fixed(w.m, pattern = "_", 2)[2],
                           condition=w.condition,
                           task=w.task,
                           repetition=w.rep,
                           dcm=w.dcm,
                           F=w.F,
                           time=w.tt,
                           iter=w.iter,
                           matrix='pH',
                           t( c( t(c(DCM@Pp@decay)),t(c(DCM@Pp@transit)),DCM@Pp@epsilon ) )  ),
               file=w.fn,
               sep=',',
               col.names=FALSE,
               row.names=FALSE,
               append=TRUE )

}






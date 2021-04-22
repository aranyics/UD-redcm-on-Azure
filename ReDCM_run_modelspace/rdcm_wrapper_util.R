
suppressPackageStartupMessages( require('stringr') )
#suppressPackageStartupMessages( require('oro.nifti') )

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

is.abs.path <- function(path) {
  # for linux
  if ( substr( path, 1, 1 ) == '/' ) {
    return( TRUE )
  }
  else {
    return( FALSE )
  }
  # TODO: for windows
  
}

colnames.cross <- function( names, sep='.' )
{
  k = length(names)
  cnames = character()
  for ( i in 1:k )
  {
    for ( j in 1:k )
    {
      n = paste( names[i], names[j], sep=sep )
      
      cnames = c( cnames, n )
    }
  }
  return(cnames)
}


maskeach.mx = function( v, u, l, miss=NULL )
{
  
  k = NULL
  j = 1:u
  if ( !is.null(miss) )
  {
    for( i in miss )
    {
      j = j[-which(j==i)]
    }
  }
  
  for ( i in j-1 )
  {
    k = c( k, v + i*l )
  }
  
  return( k )
  
}


rdcm_dec2binary = function( x )
{
  
  a = NULL
  
  while ( x>0 )
  {
    a = c( x %% 2, a); x = floor( x / 2 );
  }
  
  return( a )
  
}

rdcm_binary2matrix = function( x )
{
  
  l = length(x)
  M = diag(sqrt(l))
  x = as.logical(x)*1
  x = rev( x )
  x = x[-which(M==1)]
  k=0
  
  for ( i in 1:length(x) )
  {
    k = k + x[i] * 2^(i-1)
  }
  
  return( k )
  
}

rdcm_matrix2binary = function( x, m, nodiag=FALSE )
{
  
  l = m*m - m
  A = array( 0, l )
  M = diag(m)
  
  for ( i in 1:l )
  {
    A[i] = x %% 2; x = floor( x / 2 );
  }
  A = rev(A)
  
  M[which(M!=1)] = A
  
  if ( nodiag )
  {
    M = M - diag(m)
  }
  
  return( M )
  
}

rdcm_isolated_matrices = function( m, direct, connected=T )
{

  l = m*m - m
  interval = 0:2^l
  
#   if ( connected )
#     return( which( aaply(interval, 1, function(x){ !any(rowSums(rdcm_matrix2binary(x,m,T)) == 0 ) }) ) - 1 )
#   else
#     return( which( aaply(interval, 1, function(x){ any(rowSums(rdcm_matrix2binary(x,m,T)) == 0 ) }) ) - 1 )

  if ( connected )
    return( which( aaply(interval, 1, function(x){ !any(colSums(rdcm_matrix2binary(x,m,T))[direct] == 0 ) }) ) - 1 )
  else
    return( which( aaply(interval, 1, function(x){ any(colSums(rdcm_matrix2binary(x,m,T))[direct] == 0 ) }) ) - 1 )
  
}

rdcm_hist_highlight <- function(x, value, col.value, col=NA, ...){
  hst <- hist(x, ...)
  idx <- findInterval(value, hst$breaks, all.inside=TRUE)
  idx[which(idx == 1)] = -Inf
  cols <- rep(col, length(hst$counts))
  cols[idx] <- col.value
  hist(x, col=cols, ...)
  text(hst$mids[idx],max(hst$counts),labels=1:length(value), adj=c(0.5, -0.5), cex=0.7)
  text(hst$mids[idx],0,labels=value, adj=c(0.5, 1.2), cex=0.7)
}


rdcm_get_missing_models = function( n, u, space.dir, missB=NULL )
{
  
  #source( '/home/aranyics/R/rdcm/R/rdcm_estimate.r' )
  #source( '/home/aranyics/R/rdcm/R/rdcm_utils.r' )
  options('scipen'=10)
  
  a = diag( n )
  diagIdx = which( a == 1 )
  b = rep(0, n*n*u)
  AB = c(a, b)
  k = 1
  missing = NULL
  #df = data.frame( ID=0, A=0, B=0, B1=0, B2=0, B3=0, Fe=0, stringsAsFactors=FALSE )
  df = data.frame( ID=0, A=0, B1=0, B2=0, B3=0, Fe=0, stringsAsFactors=FALSE )
  combin = 2^(n*n-n)
  img = array(-6000, dim=c(combin, combin, combin, combin))
  
  
  for ( i in 0:(2^(n*n-n)-1) )
  {
    a = rdcm_matrix2binary(i, n)
    
    B.en = c(maskeach.mx( which( (a - diag(n)) != 0), u, n*n, missB))
    for ( j in 0:(2^(length(B.en))-1) )
    {
      b = rep(0, n*n*u)
      tmp = rdcm_dec2binary( j )
      b[B.en] = c( rep(0, length(B.en) - length(tmp)), tmp )
      
      
      if ( k %% 1000 == 0 )
        sub.dir = str_c( space.dir, "/dcm/out_", floor((k-1)/1000)*1000+1, sep='' )
      else
        sub.dir = str_c( space.dir, "/dcm/out_", floor(k/1000)*1000+1, sep='' )
      
        #if ( any(list.files(sub.dir, pattern=as.character(k), recursive=TRUE) == str_c('DCM_', k, sep='')) )
        if ( file.exists( str_c(sub.dir, '/DCM_', k, sep='') ) )
            {
#               #if ( k < 185000 )
#               #{
#                 load( str_c(sub.dir, "/DCM_", k, sep='') )
#                 if ( any(as.logical( diag(DCMe@M[[1]]@pC)[1:(n*n)] )*1 != c(a)) || any(as.logical( diag(DCMe@M[[1]]@pC)[(n*n+1):(n*n*(u+1))] )*1 != b) )
#                 {
#                  missing = c(missing, k)
#                  write.table(k, str_c(space.dir, '/missing.csv'), row.names=FALSE, col.names=FALSE, append=TRUE)
#                  write.table(str_c(a, sep='',collapse=''), str_c(space.dir, '/missingA.csv'), row.names=FALSE, col.names=FALSE, append=TRUE)
#                  write.table(str_c(b, sep='',collapse=''), str_c(space.dir, '/missingB.csv'), row.names=FALSE, col.names=FALSE, append=TRUE)
#                 }
#               #}
            }
            else
            {
              missing = c(missing, k)
              write.table(k, str_c(space.dir, '/missing.csv'), row.names=FALSE, col.names=FALSE, append=TRUE)
              write.table(str_c(a, sep='',collapse=''), str_c(space.dir, '/missingA.csv'), row.names=FALSE, col.names=FALSE, append=TRUE)
              write.table(str_c(b, sep='',collapse=''), str_c(space.dir, '/missingB.csv'), row.names=FALSE, col.names=FALSE, append=TRUE)
            }
      
#       #sub.dir = str_c( space.dir, "/ges/AB_bw" )
#       load( str_c(sub.dir, "/DCM_", k, sep='') )
#       #DCMe = DCMe.traj
#       Ai = rdcm_binary2matrix( as.logical( diag(DCMe@M[[1]]@pC)[1:(n*n)] )*1 )
#       B1i = rdcm_binary2matrix( as.logical( diag(DCMe@M[[1]]@pC)[(n*n+1):(n*n*(2))] )*1 )
#       B2i = rdcm_binary2matrix( as.logical( diag(DCMe@M[[1]]@pC)[(n*n*2+1):(n*n*(3))] )*1 )
#       B3i = rdcm_binary2matrix( as.logical( diag(DCMe@M[[1]]@pC)[(n*n*3+1):(n*n*(4))] )*1 )
#       df[k, ] = c( k,
#                    Ai,
#                    B1i,
#                    B2i,
#                    B3i,
#                    DCMe@Fe )
#       img[B1i+1, B2i+1, B3i+1, Ai+1] = DCMe@Fe
      #write.csv( df, str_c( space.dir, '/ges/AB_bw/ges_AB_bw.csv' ), row.names=FALSE )
      
      if ( k %% 10000 == 0 )
      {
        cat(k, ' ')
        #return( df )
      }
      k = k+1
    }
    
  }
  
#   dir.create( str_c( space.dir, '/fullspace' ), showWarnings=FALSE )
#   write.csv( df, str_c( space.dir, '/fullspace/fullspace_3node.csv' ), row.names=FALSE )
#   writeNIfTI( nifti(img, datatype=16), str_c( space.dir, '/fullspace/fullspace_3node' ) )
  
  return( c(k, missing) )
  
}

rdcm_get_missing_models_v2 = function( n, u, space.dir, missB=NULL )
{
  
  #source( '/home/aranyics/R/rdcm/R/rdcm_estimate.r' )
  #source( '/home/aranyics/R/rdcm/R/rdcm_utils.r' )
  options('scipen'=10)
  
  a = diag( n )
  diagIdx = which( a == 1 )
  b = rep(0, n*n*u)
  AB = c(a, b)
  k = 1
  missing = NULL
  #df = data.frame( ID=0, A=0, B=0, B1=0, B2=0, B3=0, Fe=0, stringsAsFactors=FALSE )
  df = data.frame( ID=0, A=0, B1=0, B2=0, B3=0, Fe=0, stringsAsFactors=FALSE )
  combin = 2^(n*n-n)
  img = array(-6000, dim=c(combin, combin, combin, combin))
  
  
  for ( i in 0:(2^(n*n-n)-1) )
  {
    a = rdcm_matrix2binary(i, n)
    
    B.en = c(maskeach.mx( which( (a - diag(n)) != 0), u, n*n, missB))
    for ( j in 0:(2^(length(B.en))-1) )
    {
      b = rep(0, n*n*u)
      tmp = rdcm_dec2binary( j )
      b[B.en] = c( rep(0, length(B.en) - length(tmp)), tmp )
      
      
      if ( k %% 1000 == 0 )
        sub.dir = str_c( space.dir, "/dcm/out_", floor((k-1)/1000)*1000+1, sep='' )
      else
        sub.dir = str_c( space.dir, "/dcm/out_", floor(k/1000)*1000+1, sep='' )
      
      #       if ( any(list.files(sub.dir, pattern=as.character(k), recursive=TRUE) == str_c('DCM_', k, sep='')) )
      #       {
      #         #if ( k < 185000 )
      #         #{
      #           load( str_c(sub.dir, "/DCM_", k, sep='') )
      #           if ( any(as.logical( diag(DCMe@M[[1]]@pC)[1:(n*n)] )*1 != c(a)) || any(as.logical( diag(DCMe@M[[1]]@pC)[(n*n+1):(n*n*(u+1))] )*1 != b) )
      #           {
      #            missing = c(missing, k)
      #            write.table(k, str_c(space.dir, '/missing.csv'), row.names=FALSE, col.names=FALSE, append=TRUE)
      #            write.table(str_c(a, sep='',collapse=''), str_c(space.dir, '/missingA.csv'), row.names=FALSE, col.names=FALSE, append=TRUE)
      #            write.table(str_c(b, sep='',collapse=''), str_c(space.dir, '/missingB.csv'), row.names=FALSE, col.names=FALSE, append=TRUE)
      #           }
      #         #}
      #       }
      #       else
      #       {
      #         missing = c(missing, k)
      #         #write.table(k, str_c(space.dir, '/missing.csv'), row.names=FALSE, col.names=FALSE, append=TRUE)
      #         #write.table(str_c(a, sep='',collapse=''), str_c(space.dir, '/missingA.csv'), row.names=FALSE, col.names=FALSE, append=TRUE)
      #         #write.table(str_c(b, sep='',collapse=''), str_c(space.dir, '/missingB.csv'), row.names=FALSE, col.names=FALSE, append=TRUE)
      #       }
      
    #sub.dir = str_c( space.dir, "/ges/AB_bw" )
      load( str_c(sub.dir, "/DCM_", k, sep='') )
    #DCMe = DCMe.traj
      Ai = rdcm_binary2matrix( as.logical( diag(DCMe@M[[1]]@pC)[1:(n*n)] )*1 )
      B1i = rdcm_binary2matrix( as.logical( diag(DCMe@M[[1]]@pC)[(n*n+1):(n*n*(2))] )*1 )
      B2i = rdcm_binary2matrix( as.logical( diag(DCMe@M[[1]]@pC)[(n*n*2+1):(n*n*(3))] )*1 )
      B3i = rdcm_binary2matrix( as.logical( diag(DCMe@M[[1]]@pC)[(n*n*3+1):(n*n*(4))] )*1 )
      df[k, ] = c( k,
                   Ai,
                   B1i,
                   B2i,
                   B3i,
                   DCMe@Fe )
      img[B1i+1, B2i+1, B3i+1, Ai+1] = DCMe@Fe
    #write.csv( df, str_c( space.dir, '/ges/AB_bw/ges_AB_bw.csv' ), row.names=FALSE )
      
      if ( k %% 10000 == 0 )
      {
        cat(k, ' ')
        #return( df )
      }
      k = k+1
    }
    
  }
  
  dir.create( str_c( space.dir, '/fullspace' ), showWarnings=FALSE )
  write.csv( df, str_c( space.dir, '/fullspace/fullspace_3node_C2.csv' ), row.names=FALSE )
  writeNIfTI( nifti(img, datatype=16), str_c( space.dir, '/fullspace/fullspace_3node_C2' ) )
  
  #return( c(k, missing) )
  
}

rdcm_create_bitmap = function( n, out.dir=NULL )
{
  
  #dir.create( out.dir )
  
  #df = data.frame( k = 0, b = str_c( rep(0, n), sep='',collapse=''), stringsAsFactors=FALSE )
  
  for ( k in 0:(2^n - 1) )
  {
    b = rdcm_dec2binary( k )
    write.table( t(c(k, rep(0, n-length(b)), b )), str_c(out.dir, 'test.csv'), sep=',', row.names=FALSE, col.names=FALSE, append=TRUE)
    #df[k+1, ] = c( k, str_c(c( rep(0, n-length(b)), b ), sep='',collapse='') )
  }
  
  #return( df )
  
}

#Read in RCT3-style files



#String to array of strings
sepfortranstring = function(basestring)
{
  #Treats multiple delimiters as a single delimiter
  #Strings must not have spaces in

  if(!is.character(basestring))
  {
    warning("basestring is not a character variable, converting with 'as.character'")
    basestring = as.character(basestring)
  } 
  #replace tabs with spaces
  basestring = paste (unlist(strsplit(basestring, "\t")), collapse=" ")  
  #replace commmas with spaces
  basestring = paste (unlist(strsplit(basestring, ",")), collapse=" ")
  
  #remove extraneous quotes from within strings
  basestring = paste (unlist(strsplit(basestring, "\"")), collapse="")
  retval = unlist(strsplit(basestring, " "))
  retval[retval!=""]
  
}








readRCT3 = function(filename,old=TRUE)
{
  ow <- options()$warn
  options(warn = -1)   ##Not while debugging!
  on.exit(options(warn = ow))

  head1 = readLines(filename,n=1)
  head2 = sepfortranstring(scan(filename,what=character(1),skip=1,n=1,sep="\n",quiet=TRUE))
  head2 = as.numeric(head2[1:3])
  if (is.na(head2[1]) | is.na(head2[3]) | is.na(head2[3]))
  {    
    cat("cannot read file information from header:", #filename, 
      "", sep="\n")
    return(NULL)
  } 
  nsurveys = head2[1]
  nyears = head2[2]
  VPAn = head2[3]
  filecontents = data.frame(matrix(0,nyears,nsurveys+2))
  if (old)
  { #Old file format
    for (year in 1:nyears) filecontents[year,1:(nsurveys+2)] = as.numeric(sepfortranstring(scan(filename,what=character(1),skip=year+1,n=1,sep="\n",quiet=TRUE)))[1:(nsurveys+2)]       
    colnames(filecontents)[1:2] = c("Year","VPA")
    for (survey in 1:nsurveys) colnames(filecontents)[survey+2] = sepfortranstring(scan(filename,what=character(1),skip=nyears+survey+2,n=1,sep="\n",quiet=TRUE))[1]       
    if (is.na(sum(filecontents))) 
    { 
      cat("File can't be read as an old style file correctly\n")
      return(NULL)
    }
  } else { #New file format
    for (year in 1:nyears) filecontents[year,1:(nsurveys+2)] = as.numeric(sepfortranstring(scan(filename,what=character(1),skip=year+1,n=1,sep="\n",quiet=TRUE)))[1:(nsurveys+2)]       
    colnames(filecontents)[1:2] = c("Year","VPA")
    colnames(filecontents)[3:(nsurveys+2)] = sepfortranstring(scan(filename,what=character(1),skip=2,n=1,sep="\n",quiet=TRUE))[3:(nsurveys+2)]      
    if (is.na(sum(filecontents)))
    {
      cat("File can't be read as an new style file correctly\n")
      return(NULL)
    }
  }
  
  ##What if VPAn != 2
  if (VPAn > 2)  filecontents = filecontents[,c(1,VPAn,2:(VPAn-1),(VPAn+1):(nsurveys+2))]
  if (VPAn == 1) filecontents = filecontents[,c(2,1,3:(nsurveys+2))]  
  print(paste("File details:",head1))
  list(filecontents=filecontents, head=head1)

}


# IRISH SEA PLAICE CONFIGURATION
# 'C'  Regression type
# 'N'  Tapered weighting
# 0    Taper power
# 0    Taper range
# 'Y'  Shrink to mean
# 'N'  Exclude surveys with SE greater than mean
# 5    Minimum number of regression points
# 0.2  Minimum S.E.

readConfigFile <- function(filename,  old, configfile)
{
  filetext = readLines(configfile)
  regtype =     strsplit(temp,"")[[2]][2]
  taperWeight = strsplit(temp,"")[[3]][2]
  taperPower =  as.numeric(strsplit(temp," ")[[4]][1])
  taperRange =  as.numeric(strsplit(temp," ")[[5]][1])
  shrink =      ifelse(strsplit(temp,"")[[6]][2] %in% c("y","Y"), TRUE, FALSE)
  exclude  =    strsplit(temp,"")[[7]][2]
  minPoints =   as.numeric(strsplit(temp," ")[[8]][1])
  minSE =       as.numeric(strsplit(temp," ")[[9]][1])
  print(paste("Config file:",filetext[1]))
  if (taperWeight %in% c('n','N'))
  {
    RCT3(filename, old=old, logged = NA, minpoints = minPoints, regression=regtype, shrink=shrink, taper=0, tapert=9999, minse=minSE, configfile = NA)
  } else {
    RCT3(filename, old=old, logged = NA, minpoints = minPoints, regression=regtype, shrink=shrink, taper=taperpower, tapert=taperRange, minse=minSE, configfile = NA)  
  }
}

#Like sprint except first argument sets maximum string length, longer strings are 
#replaced with "*"s as in Fortran output. NAs => "  "
sprintf2 <- function(n, ...)
{
   temp = sprintf(...)
   #if (nchar(temp) > n) temp = nstr("*",n)
   #if (temp == paste(nstr(" ",nchar(temp)-2),"NA",sep="")) temp = nstr(" ",n)
   temp
}

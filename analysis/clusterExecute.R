



ghpc <<- "aw30w@ghpcc06.umassrc.org"
proj.hpc <<- "/project/umw_zhiping_weng/wespisea/"
home.hpc <<- "/home/aw30w"
zlab.scripts <<- "/home/wespisea/log/rscripts/"
ghpc.scripts <<- "/home/aw30w/log/rscripts/"


shCheckFile <- function(f){
  paste0("[[ -s ",f," ]] && echo TRUE || echo FALSE")
}


hpc.file.exists <- function(f,...){
  as.logical(hpc.system(shCheckFile(c(f,...))))
}

scpFile <- function(file.local, dir.remote){
  if(!file.exists(file.local)){
    stop("cannot transfer file that doesn't exist")
  }
  
  cmd <- paste0("scp ",file.local," " ,ghpc,":",dir.remote)
  system(cmd)
}

hpc.system.helper <- function(cmd){
  r.cmd <- pipe(paste("ssh",ghpc,"\"",cmd,"\""))
  str <- readLines(r.cmd)
  flush(r.cmd)
  close(r.cmd)
  rm(r.cmd)
  str
}
hpc.system <- function(cmd){
  sapply(cmd,hpc.system.helper)  
}

# example with $ : hpc.system("echo \\$HOME") -> /home/aw30w

hpc.system.nowait <- function(cmd){
  r.cmd <- paste("ssh",ghpc,"\"",cmd,"\"")
  system(r.cmd,wait=FALSE)
}

hpc.file.exists.old <- function(file){
  cmd <- paste("file ",file)
  result <- hpc.system(cmd)
  g.result <- grep(x=result, pattern="No such file or directory")
  if(any(g.result)){
    return(FALSE)
  } else {
    return(TRUE)
  }
}

hpc.file.exists <- function(f,...){
  as.logical(hpc.system(shCheckFile(c(f,...))))
}


hpc.file.md5sum.helper <- function(file){
  if(!hpc.file.exists(file)){
    return(NA)
  }
  cmd <- paste("md5sum",file,"|cut -d ' ' -f 1")
  result <- hpc.system(cmd)
}

hpc.file.md5sum <- function(file){
  sapply(file, function(f)hpc.file.md5sum.helper(f))
}



# url.vec <- c("cnn.com", "facebook.com", "bbc.co.uk")
# target.vec <- paste0("~/sandbox/wget/",c("cnn", "facebook", "bbc"))
hpc.download <- function(url.vec,target.vec){
  if(missing(url.vec) || missing(target.vec) || (!identical(length(target.vec),length(url.vec)))){
    stop("must have url.vec and target.vec of same length")
  }
  cmd <- paste(paste0(paste0("wget ", url.vec, " -O ", target.vec),collapse= " & \n"), "&")
  hpc.system.nowait(cmd)
}


# url.vec <- c("cnn.com", "facebook.com", "bbc.co.uk")
# target.vec <- paste0("~/sandbox/wget/",c("cnn", "facebook", "bbc"))
hpc.download.seq.cont <- function(url.vec,target.vec){
  if(missing(url.vec) || missing(target.vec) || (!identical(length(target.vec),length(url.vec)))){
    stop("must have url.vec and target.vec of same length")
  }
  cmd <- paste0(paste0("wget --continue ", url.vec, " -O ", target.vec),collapse= " ; \n")
  hpc.system.nowait(cmd)
}





createBsubCmd <- function(cmd,cores=4,memory=4,wait=60,queue="short",abbrv="testScript"){

memoryGB = memory * 1024  
abbr <- paste0(abbrv,"-",format(Sys.time(),"%m%d-%H%M%S"))
local.file <- file.path(zlab.scripts,paste0(abbr,".sh"))
if (file.exists(local.file)){
  abbr <- paste0(abbrv,"-",format(Sys.time(),"%m%d-%H%M%S"),"r",gsub(gettextf("%3.f",runif(1)*1000), pattern=" ",replacement=""))
  local.file <- file.path(zlab.scripts,paste0(abbr,".sh"))
}
print(local.file)
remote.file <- file.path(ghpc.scripts,paste0(abbr,".sh"))
remote.err.file <- file.path(ghpc.scripts,paste0(abbr,".err"))
remote.out.file <- file.path(ghpc.scripts,paste0(abbr,".out"))

scriptstr <- gettextf("#!/bin/sh
#BSUB -L /bin/bash
#BSUB -n %i 
#BSUB -R \"rusage[mem=%i] span[hosts=1]\"
#BSUB -q %s
#BSUB -W %i
#BSUB -e %s 
#BSUB -o %s

%s 
",cores,memoryGB,queue,wait,remote.err.file,remote.out.file,cmd)

cat(scriptstr,file=local.file)
list(local=local.file,remote=remote.file,error=remote.err.file,out=remote.out.file)

}


countMappedReads <- function(bamFile){
hpc.system(paste("bash;source ~/.bashrc;/home/aw30w/bin/countMappedReads ",bamFile))
}

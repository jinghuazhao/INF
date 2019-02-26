# Manhattan plot script
# Basic script written by Joshua C Randall & Reedik Magi, further adjusted by BP

# Read in arguments from the command line

rm(list=ls())
for (e in commandArgs(trailingOnly=TRUE))
{
  ta = strsplit(e,"=",fixed=TRUE)
  if(!is.null(ta[[1]][2]))
  {
    assign(ta[[1]][1],ta[[1]][2])
  } else {
    assign(ta[[1]][1],TRUE)
  }
}
if(!exists("input"))
{
  input <- paste("snptest.out")
}
if(!exists("output")) {
  output <- paste(input,".manh.png",sep="")
}

### Defining the output image type
png(output,width=3600, height=2400, pointsize = 12, res=300)

# Define significance level
p_sign<-5e-8

### Reading in GWAS datafile and known variants

#stub<-"GDNF.unrelated_4994_pihat_0.1875_autosomal_typed_only.snptest.plot-data.out"
#data<-fread(paste(stub,sep=""),na.strings = c("-9","-1","NA","-nan")) # Ask Praveen to help installing the data.table package, doesn't work on cardio for this package
data<-read.table(gzfile(input),na.strings = c("-9","-1","NA","-nan"),stringsAsFactors=FALSE,header=TRUE)

# Define plot title
# Remove the "/", and get the first string when do string split on "_", which is the pheno name (${PHENOTYPE_NAME}_interval_olink_subset_unrelated_4994_pihat_0.1875_autosomal_typed_only_)
input_file_path<-strsplit(input,"/",fixed=TRUE)
plot_title<-(strsplit(input_file_path[[1]][(length(input_file_path[[1]]))],"_",fixed=TRUE))[[1]][1]

### Work the data a bit
colnames(data)<-c("chromome","position","p_value")
tmp_store<-data
data<-data[complete.cases(data),]  

obspval <- (data$p_value)
chr <- (data$chromome)
pos <- (data$position)


### Define plot dimensions
# How many chromosomes ?
max_nchr<-length(unique(chr))


# Defining the maximum p-value to set the Y-axis size / Setting -10 Log P-values Y-axis max above 250 to 250
obspval.sane<-obspval[!is.na(obspval) & (obspval > 0)]
obspval.clean<-ifelse(obspval.sane < 300, obspval.sane, 300) 
obsmax.trunc <- trunc(max(-log10(obspval.clean)))
obsmax.true<-max(-log10(obspval.clean))
# Nice truncations
obsmax.ylim<-ifelse(obsmax.trunc < 10,10, 
ifelse(obsmax.trunc < 15, 15,
ifelse(obsmax.trunc < 15, 15,
(((trunc(obsmax.trunc/10)+1)*10)))))


# Sorting the plotting data from the GWAS datafile
sort.ind <- order(chr, pos) 
chr <- chr[sort.ind]
pos <- pos[sort.ind]
obspval <- obspval[sort.ind]


### Setting the plotting data
# Set indices
x <- 1:max_nchr
x2<- 1:max_nchr

# Re-assign chromosome numbers not based on real chr nr, but based on order, so that there won't be gaps
max_nchr<-length(unique(chr))
chr_index<-unique(chr)
chr_new_index<-1:max_nchr

for (i in 1:max_nchr)
{
chr[which(chr==chr_index[i])]<-rep(chr_new_index[i],length(chr[which(chr==chr_index[i])]) )
}

# Define the actual plotting coordinates
for (i in 1:max_nchr)
{
	 curchr=which(chr==i)
	 x[i] <- trunc((max(na.omit(pos[curchr])))/100) +100000
	 x2[i] <- trunc((min(na.omit(pos[curchr])))/100) -100000
}

x[1]=x[1]-x2[1]
x2[1]=0-x2[1]

for (i in 2:(max_nchr+1))
{
	x[i] <- x[i-1]-x2[i]+x[i]
	x2[i] <- x[i-1]-x2[i]
}

locX = trunc(pos/100) + x2[chr]
locY = -log10(obspval)




# Set the default colors of the datapoints, grey and light grey, alternating between odd and even chromosome numbers
col1=rgb(176,176,176,maxColorValue=255)
col2=rgb(158,158,158,maxColorValue=255)
curcol <- ifelse (chr%%2==0, col1, col2)

### Start plotting the data
plot(locX,locY,pch=20,col=curcol,axes=F,ylab="",xlab="",bty="n",ylim=c(0,obsmax.ylim),cex=0.8,main=plot_title)

# Setting horizontal reference lines, uncomment if they need to be plotted, genomide=height of reference line
p_sign_thresh_plot=-log10(p_sign)
lines(c(0,max(locX)),c(p_sign_thresh_plot,p_sign_thresh_plot),col="dodgerblue4", lty=2,lwd=1)

### If there are hits that are below <p_sign>, highlight these peaks
if (obsmax.true>=-log10(p_sign))
{

# Define the location of significant loci (smaller than -log10(<p_sign>)) that need to be higlighted
highlight_flank<-10000

hitsposition=  locX[locY >= (-log10(p_sign))]

## tophits location +/- <highlight_flank> for highlighting peaks
locXmin=hitsposition-highlight_flank
locXmax=hitsposition+highlight_flank

# Define colors for seperate hits, they will be colored blue / read alternatively
color_defs=c("navy","firebrick3")
color_vec<-rep(color_defs,length(hitsposition))
hits_color_vector<-NULL
hits_position_min<-hitsposition-highlight_flank
current_hits_position_min=min(hits_position_min)
current_color<-1

# If there are more than 1 seperate hits, color these blue and red alternatively
for (i in 1:length(hits_position_min)) 
{

if (hits_position_min[i]-current_hits_position_min<=(2*highlight_flank))
{
current_hitsposition_min<-current_hits_position_min
hits_color_vector<-c(hits_color_vector,color_vec[current_color])
}
else
{
current_hits_position_min<-hits_position_min[i]
current_color<-current_color+1
hits_color_vector<-c(hits_color_vector,color_vec[current_color])
}
}

# Plot the hit points blue an dred
for (i in 1:length(hitsposition)) 
{
points(locX[which(locX > locXmin[i] & locX < locXmax[i])], locY[which(locX > locXmin[i] & locX < locXmax[i])], col=hits_color_vector[i], pch=20, cex=0.8)
}
}

# Draw the y-axis
axis(2,las=1,pos=0)

# Make (non-existing) x-axis text and labels
default.full.x.axis.txt<-c(1:22,"X","Y XY M","","")
plot.x.axis.txt<-default.full.x.axis.txt[chr_index]

for (i in 1:max_nchr)
{
	labpos = (x[i] + x2[i]) / 2
	mtext(plot.x.axis.txt[i],1,at=labpos,cex=1,line=0)
}

# Set the horizontal axis text - "Chromosome" and each chromosome number for which the data is plotted
mtext("Chromosome",1,at=x[max_nchr]/2,cex=1,line=2)
mtext("-log10 p-value",2,line=1)

# Finally, make a legend
# legend('topright', as.character(c("Novel loci exome-wide significant       ","Known loci exome-wide significant       ")), pch=c(19,19), col=c("firebrick3","navy"), bty='n', cex=1)

# Stop imaging device
dev.off()

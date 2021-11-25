######################################################################################################################################################
###                                                                                                                                                ###
###     turboqq : An R script to produce QQ plots for impatient people. Concentration bands code from an orginal script by Mike Weale and          ###
###              Tom Price, King's College London : https://sites.google.com/site/mikeweale/software/qqplots/qq_plot_v7.r                          ###
###              Concentration bands are plotted using the pointwise method of Quesenberry & Hale (1980) J. Statist. Comput. Simul. 11:41-53       ###
###              The method proceeds from noting that the kth order statistic from a sample of n i.i.d. U(0,1) statistics has a Beta(k,n+1-k)      ###
###              distribution.                                                                                                                     ###
###                                                                                                                                                ###
###     Version : 0.1.0  : BETA VERSION WITHOUT WARRANTIES !                                                                                       ###
###                                                                                                                                                ###
######################################################################################################################################################

#=====================================================================
# Changelog 
#=====================================================================
# 2017/XX/XX : ----
#=====================================================================

#---------------------------------------------------------------------
# Running the script from command line
#---------------------------------------------------------------------
# R --slave --vanilla --args \
# input_data_path=${PWD}/my_input_assoc_data_file \
# output_data_rootname=${PWD}/my_qq_plot \
# plot_title="my plot title" < turboqq.r
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Input association data path / input_data_path
#---------------------------------------------------------------------
# Define path of the input association data
# The input data needs to be a file that has :
# 1. Spaces as field separators
# 2. One header line
# 3. Option I  (no extreme p-values present): 3 columns, being 
#               chromsome, position, pvalue - in this order, 
#               column names are not important
#    Option II (extreme p-values present): 5 columns, being
#               chromsome, position, pvalue, beta, se 
#               - in this order, column names are not important
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Output data rootname / output_data_rootname
#---------------------------------------------------------------------
# Define root name of the plot output file
# Ex. : "my_qq_plot" will result 
# in an output file named "my_qq_plot.png"
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Plot title / plot_title
#---------------------------------------------------------------------
# Define plot title which will be displayed on top of the plot
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# To do list
#---------------------------------------------------------------------
# 1. ---
#---------------------------------------------------------------------
# Read in arguments from the command line

rm(list=ls())

## Verbose printing status bars
fat_status_bar<-"============================================================================================================"
skinny_status_bar<-"------------------------------------------------------------------------------------------------------------"

## Attention
print("")
print(fat_status_bar)
print(" 1. Read in arguments from the command line")
print(fat_status_bar)
print("")


for (arg in commandArgs(trailingOnly=TRUE))
{
  ta = strsplit(arg,"=",fixed=TRUE)
  if(!is.na(ta[[1]][2]))
  {
    assign(ta[[1]][1],ta[[1]][2])
  } else {
    stop("Not all arguments are given")
  }
}

print(paste0("  Data file path : ",input_data_path))

## Assign variable classes
input_data_path <- as.character(input_data_path)
output_data_rootname <- as.character(output_data_rootname)
plot_title <- as.character(plot_title)


#---------------------------------------------------------------------
# Reading in association plot data with scan
#---------------------------------------------------------------------

## Attention
print("")
print(fat_status_bar)
print(" 2. Reading in association plot data with scan")
print(fat_status_bar)
print("")




initial_data_dims<-dim(as.data.frame(read.table(gzfile(input_data_path), header=TRUE, stringsAsFactors=FALSE, nrows=10)))[2]

if (initial_data_dims==3) {
initial_data <- data.frame(scan(input_data_path,
                        what = list(chromosome = 0, position = 0,pvalue= 0),
                        skip=1,
                        sep=" ",
                        quiet=TRUE))

initial_data_contains_beta_se<-FALSE

} else if (initial_data_dims==5) {
initial_data <- data.frame(scan(input_data_path,
                        what = list(chromosome = 0, position = 0,pvalue= 0, beta=0, se=0),
                        skip=1,
                        sep=" ",
                        quiet=TRUE))
                        
initial_data_contains_beta_se<-TRUE

} else {

  stop("Input data does not have expected dimensions")
  }

## Attention
print("")
print(fat_status_bar)
print(" 3. Calculate log P values")
print(fat_status_bar)
print("")

  
## Check if p-values are already logged
if (length(which(initial_data$pvalue>1))>0) {
    
    initial_data$log_pvalue<-initial_data$pvalue
    ## Get only the complete data
    initial_data<-initial_data[complete.cases(initial_data),]

    # Remove the original pvalues
    initial_data$pvalue<-NULL
} else {

    ## Calculate the -log10 p-value for the input data 
    initial_data$log_pvalue<--log10(initial_data$pvalue)

    ## If beta/SE are provided, and pvalues are missing (because they are extreme), log10 P recalculate from beta/SE
    missing_pvalues_index<-which((is.na(initial_data$log_pvalue) | initial_data$log_pvalue==0))
    
    if (initial_data_contains_beta_se & (length(missing_pvalues_index)>0)) {
        # Calculate expected p-values for missing data
        missing_pvalues<-(-log(2, base=10)-pnorm(-abs(initial_data$beta[missing_pvalues_index]/initial_data$se[missing_pvalues_index]), log=T)/log(10))
        # Only replace if indeed they were below the smallest non-zero normalized floating-point number
        initial_data[missing_pvalues_index,c("log_pvalue")]<-ifelse(missing_pvalues > -log10(.Machine$double.xmin),missing_pvalues,NA)
    }

    ## Get only the complete data
    initial_data<-initial_data[complete.cases(initial_data),]

    # Remove the original pvalues
    initial_data$pvalue<-NULL
}

#-----------------------------------------------------------------------------------------
# Make settings
#-----------------------------------------------------------------------------------------

# Plot resolution
plot_resolution<-1800

## For calculating concentration bands
alpha=0.05 # for calculating concentration bands
df=1
one.sided=FALSE

#-----------------------------------------------------------------------------------------
# Calculate plotting parameters
#-----------------------------------------------------------------------------------------

## Attention
print("")
print(fat_status_bar)
print(" 4. Calculate plotting points")
print(fat_status_bar)
print("")


## Sort pvals and calculate expected pvals
log_obspval_sorted <- sort(initial_data$log_pvalue,decreasing=TRUE)
exppval <- c(1:length(log_obspval_sorted))
log_exppval <- -(log10( (exppval-0.5)/length(exppval)))

## Get the maxima
ymax=max(na.omit(log_obspval_sorted))
xmax=max(na.omit(log_exppval))
    
# Calculate the lambda  
lambda <- median(qchisq(na.omit((log_obspval_sorted*-1*log(10))), df=1, lower.tail=FALSE, log.p = TRUE)) / qchisq(0.5, 1)

## Set vertical resolution
# We will maximally allow a fixed number of points to be plotted vertically, 
# here choosing 800 as a 'pixel' unit on a normal standard R plot.
vertical_resolution<-plot_resolution
horizontal_resolution<-plot_resolution

## Now we will scale the 800-point resolution for the p-values
obs_log_pvalue_break_size<-ymax/vertical_resolution
exp_log_pvalue_break_size<-xmax/horizontal_resolution

## Create a vector from 0 to the vertical resolution, which we will use to bin pvalues
obs_log_pvalue_scaling_vector<-seq(0,vertical_resolution,by=obs_log_pvalue_break_size)
exp_log_pvalue_scaling_vector<-seq(0,horizontal_resolution,by=exp_log_pvalue_break_size)

## Bin the pvals and get only unique pairs
obs_log_pvalue_binned<-.bincode(log_obspval_sorted, obs_log_pvalue_scaling_vector, right = TRUE, include.lowest = FALSE)*obs_log_pvalue_break_size
exp_log_pvalue_binned<-.bincode(log_exppval, exp_log_pvalue_scaling_vector, right = TRUE, include.lowest = FALSE)*exp_log_pvalue_break_size
plot_data_binned<-as.data.frame(cbind(exp_log_pvalue_binned,obs_log_pvalue_binned))
plot_data_reduced<-unique(plot_data_binned)

## Truncate the maxima and use these to make  nice ticks
pretty_ymax=trunc(ymax+1)
pretty_bigymax=trunc(ymax)
pretty_xmax=trunc(xmax+1)
if (ymax<=15) {
y4Ly = c(0:pretty_ymax)
    } else if ((ymax>15) & (ymax<=100)) { 
    y4Ly = c(seq(0,pretty_bigymax+10,10))
    } else if ((ymax>100) & (ymax<=200)) { 
    y4Ly = c(seq(0,pretty_bigymax+20,20)) 
    } else if ((ymax>200) & (ymax<=300)) { 
    y4Ly = c(seq(0,pretty_bigymax+30,30)) 
    } else if ((ymax>300) & (ymax<=400)) { 
    y4Ly = c(seq(0,pretty_bigymax+40,40)) 
    } else if ((ymax>400) & (ymax<=500)) { 
    y4Ly = c(seq(0,pretty_bigymax+50,50)) 
    } else if ((ymax>500) & (ymax<=600)) { 
    y4Ly = c(seq(0,pretty_bigymax+60,60)) 
    } else if ((ymax>600) & (ymax<=700)) { 
    y4Ly = c(seq(0,pretty_bigymax+70,70)) 
    } else if ((ymax>700) & (ymax<=800)) { 
    y4Ly = c(seq(0,pretty_bigymax+80,80)) 
    } else if ((ymax>800) & (ymax<=900)) { 
    y4Ly = c(seq(0,pretty_bigymax+90,90)) 
    } else if ((ymax>900) & (ymax<=100)) { 
    y4Ly = c(seq(0,pretty_bigymax+100,100)) 
    } else { y4Ly = c(seq(0,pretty_bigymax+200,200)) 
}

x4Lx = c(0:pretty_xmax)
    xnums = (x4Lx)                   
    ynums = (y4Ly)                       
    Lx <- parse( text=paste(x4Lx,sep="") )
    Ly <- parse( text=paste(y4Ly,sep="") )

## Use the maxima for making x and y limits for plotting function
    x.lim=xmax
    y.lim=max(y4Ly)
    
#-----------------------------------------------------------------------------------------
# Start the actual plotting
#-----------------------------------------------------------------------------------------

## Attention
print("")
print(fat_status_bar)
print(" 5. Start the actual plotting")
print(fat_status_bar)
print("")

# Start device
png(paste0(output_data_rootname,".png"),height=plot_resolution,width=plot_resolution, pointsize = 12, res=300)

# Create empty plot
# Just plots the outside box
plot(0, main = plot_title, xlab = "Expected p-value", ylab = "Observed p-value", type = "n", xlim=c(0,x.lim), ylim=c(0,y.lim),xaxt='n',yaxt='n')       #Just plots the outside box
  
# Draw axes
axis(1, at=xnums, labels=Lx )
axis(2, at=ynums, labels=Ly )

# Define and draw the concentration bands
# Note that conc band won't draw if x has too many datapoints

n <- length(log_exppval)
frac=1
starti = floor((n-1)*(1-frac)) +1                            #i for the first sorted datapoint to be plotted.
lena = n-starti+1                                            #Number of datapoints to be plotted
a2=(1:lena)                                                  #indices to be plotted
b <- n+1-a2                                                  #indices used in determining concentration band 

## Attention
print("")
print(skinny_status_bar)
print("  a. Calculate and plot concentration band points")
print(skinny_status_bar)
print("")

# Define and draw the concentration bands
# Note that conc band won't draw if x has too many datapoints
if (one.sided==FALSE) {
upper <- -log10(qbeta( 1-alpha/2, a2, b ))      #Exp. upper CL for 'a'th U(0,1) order statistic (becomes 'lower')
lower <- -log10(qbeta( alpha/2, a2, b ))        #Exp. lower CL for 'a'th U(0,1) order statistic (becomes 'upper')
} else {
upper <- rep(1,length(log_exppval))                    #Exp. upper CL for 'a'th U(0,1) order statistic (becomes 'lower')
lower <- qbeta( alpha, a2, b )          #Exp. lower CL for 'a'th U(0,1) order statistic (becomes 'upper')
}


polygon( c( log_exppval, rev(log_exppval) ), c(upper, rev(log_exppval) ), col="grey", border = NA ) #'lower' band
polygon( c( log_exppval, rev(log_exppval) ), c(lower, rev(log_exppval) ), col="grey", border = NA ) #'upper' band

# =========================== CODE FOR REDUCING AND PLOTTING CONCENTRATION BANDS - VERY SLOW STILL =========================================================== #
# # Bin the pvals for upper and lower
# upper_x_binned<-.bincode(c( log_exppval, rev(log_exppval) ),exp_log_pvalue_scaling_vector, right = TRUE, include.lowest = FALSE)*exp_log_pvalue_break_size
# upper_y_binned<-.bincode(c(upper, rev(log_exppval)),exp_log_pvalue_scaling_vector, right = TRUE, include.lowest = FALSE)*exp_log_pvalue_break_size
# upper_binned<-unique(cbind(upper_x_binned,upper_y_binned))

# lower_x_binned<-.bincode(c( log_exppval, rev(log_exppval) ),exp_log_pvalue_scaling_vector, right = TRUE, include.lowest = FALSE)*exp_log_pvalue_break_size
# lower_y_binned<-.bincode(c(lower, rev(log_exppval) ),exp_log_pvalue_scaling_vector, right = TRUE, include.lowest = FALSE)*exp_log_pvalue_break_size
# lower_binned<-unique(cbind(lower_x_binned,lower_y_binned))

# # Draw the CL bands
# polygon( upper_binned[,1], upper_binned[,2], col="grey", border = NA ) #'lower' band
# polygon( lower_binned[,1], lower_binned[,2], col="grey", border = NA ) #'upper' band
# ============================================================================================================================================================ #

# Draw the diagonal
lines(c(0,xmax),c(0,xmax),col="blue", lty=2,lwd=0.5)

# Print the lambda
# text(0.120 * xmax, 0.95 * ymax, substitute(paste(lambda[GC], "=", x), list(x=formatC(round(lambda,3),3,format="f"))), col = ifelse(lambda > 1.1, "red", "black"))
legend("topleft", legend=substitute(paste(lambda[GC], "=", x), list(x=formatC(round(lambda,3),3,format="f"))),
       text.col = ifelse(lambda > 1.1, "red", "black"), bty = "n")

## Attention
print("")
print(skinny_status_bar)
print("  b. Plot points P values")
print(skinny_status_bar)
print("")

# Finally, plot points
points(plot_data_reduced[,1],plot_data_reduced[,2], pch=19, cex=0.5, col="dodgerblue4" ) 

# Turn off the device
dev.off()

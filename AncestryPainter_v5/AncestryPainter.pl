## run AncestryPainter
## Modified by Qidi Feng
## Date: 2020-01-11 14:56
#!usr/bin/env perl -w
use strict;
use Getopt::Std;
my ($inped, $inQ, $inColor, $inPoporder, $outfile, $target_pop, $if_line, $format, %opts);
die "Usage: perl AncestryPainter.pl -i ./example/data.ind -q ./example/data.Q \nMore options:
‘-q’: (required, input file) Define the individual Q-matrix file. Check ./example/data.Q for input format.

‘-i’:  (required, input file) Define the population information of each individual. Check ./example/data.ind for input format.

‘-t’: (optional, population name) Define the only target population to be plotted as a pie chart in the center of the circle figure. The population must be included in both of the required files.

‘-p’: (optional, input files) Define the populations to be included in the figure, and the display order of the populations in the figure. Check ./example/pop.order for input format.

‘-c’: (optional, input files) Define the colors of each ancestry in the figure. Check ./example/color.file for input format.

‘-o’: (optional) Define the suffix of output files. ‘out’ by default. 

‘-l nolines’: (optional) remove the black lines between populations. Default is keeping black lines (Figure6). 

‘-f png’: (optional) output png figure instead of pdf figure. Default is pdf figure.\n" if @ARGV < 2;
getopts ("i:q:c:p:o:t:l:f:", \%opts);

#####check if color file exists##############
my (@colors);
if(defined $opts{'c'}){
    $inColor = $opts{'c'};
    open INCOLOR, "$inColor";
    while(<INCOLOR>){
        $_ =~ s/^(\s+)//g; $_ =~ s/(\s+)$//g;
        my $color = (split /\s+/, $_)[0];
        push @colors, $color;
    }
    close INCOLOR;
}
###otherwise use default color setting###
else{
    print "Using default colors\n";
    @colors = qw/RoyalBlue ForestGreen gold indianred1 olivedrab1 springgreen sienna skyblue gray tomato4 darkcyan orangered dodgerblue olivedrab4 tan midnightblue orchid3 green4 indianred3 white aliceblue antiquewhite antiquewhite1 antiquewhite2 antiquewhite3 antiquewhite4 aquamarine aquamarine1 aquamarine2 aquamarine3 aquamarine4 azure azure1 azure2 azure3 azure4 beige bisque bisque1 bisque2 bisque3 bisque4 black blanchedalmond blue blue1 blue2 blue3 blue4 blueviolet brown1 brown2 brown3 brown4 burlywood burlywood1 burlywood2 burlywood3 burlywood4 cadetblue cadetblue1 cadetblue2 cadetblue3 cadetblue4 chartreuse chartreuse1 chartreuse2 chartreuse3 chartreuse4 chocolatelchocolate1 chocolate2 chocolate3 chocolate4 coral coral1 coral2 coral3 chocolate brown/;
}
#####check if inped file exists##############
if(defined $opts{'i'}){
    $inped = $opts{'i'};
}
else{
    print "Error: No ped file indicated\n";
    last;
}
#####check if inQ file exists###############
if(defined $opts{'q'}){
    $inQ = $opts{'q'};
}
else{
    print "Error: No Q matrix file indicated\n";
    last;
}
#####check if outfile name exits############
if(defined $opts{'o'}){
    $outfile = $opts{'o'};
    print "Using $outfile as the suffix of outfiles\n";
}
else{
    print "Using 'out.' as the suffix of outfiles\n";
    $outfile = 'out';
}
#####check if target_pop exists#############
my $target_pop_signal = 0;
if(defined $opts{'t'}){
    $target_pop = $opts{'t'};
    $target_pop_signal = 1;
}
else{
    print "Warning:there is no target population indicated to be drawn in the circle\n";
}

#####save population information of each individual#############
my (%pop_ind_line);
open INPED, "$inped";
while(<INPED>){
    $_ =~ s/^(\s+)//g; $_ =~ s/(\s+)$//g;
    my ($pop, $ind) = (split /\s+/, $_)[0, 1];
    $pop_ind_line{$pop}{$ind} = $.;
}
close INPED;

if(defined $opts{'t'}){
    if(not exists $pop_ind_line{$target_pop}){
        die "Error: There is no $target_pop in the ind file\n";
    }
}

#####save the admixture proportion of each individual#############
my (%qhash);
open INQ, "$inQ";
while(<INQ>){
    $_ =~ s/^(\s+)//g; $_ =~ s/(\s+)$//g;
    $qhash{$.} = $_;
}
close INQ;

    my @pop_order; 
    my %rep_pop_K;

#####order population according to its admixture proportion#############
    my (%total_k_pop, %pop_average_k, %K_prop_pop);
    for my $pop (sort keys %pop_ind_line){
        for my $ind (sort keys %{$pop_ind_line{$pop}}){
            my @proportions = split /\s+/, $qhash{$pop_ind_line{$pop}{$ind}};
            for my $k (0..$#proportions){
                $total_k_pop{$k}{$pop} += $proportions[$k];
            }
        }
    }
    for my $k (sort keys %total_k_pop){
        for my $pop (sort keys %{$total_k_pop{$k}}){
            my @ind_no = keys %{$pop_ind_line{$pop}};
            my $ind_no = $#ind_no + 1;
            my $average = $total_k_pop{$k}{$pop}/$ind_no;
            $pop_average_k{$pop}{$average} = $k;
        }
    }
    for my $pop (sort keys %pop_average_k){
        my @values = reverse sort {$a<=>$b} keys %{$pop_average_k{$pop}};
        my $representativeK = $pop_average_k{$pop}{$values[0]};
        $rep_pop_K{$pop} = $representativeK;
        $K_prop_pop{$representativeK}{$values[0]}{$pop} = 1;
    }
##check if population order file exists#######
##order population according to our default sorting if no population order file exists#######
if(defined $opts{'p'}){
    $inPoporder = $opts{'p'};
    open INPOP, "$inPoporder";
    while(<INPOP>){
        $_ =~ s/^(\s+)//g; $_ =~ s/(\s+)$//g;
        push @pop_order, (split /\s+/, $_)[0];
    }
    close INPOP;
}
else{
    print "No population order file indicated, automatically sort populations according to ancestry proportions\n";
    for my $K (sort {$a<=>$b} keys %K_prop_pop){
        for my $value (reverse sort {$a<=>$b} keys %{$K_prop_pop{$K}}){
            for my $pop (sort keys %{$K_prop_pop{$K}{$value}}){
                push @pop_order, $pop;
            }
        }
    }
}
#### sort individuals based on their major component proportion for each population ###
print "Producing .ancestry file";
my $ancestry_file = $outfile.'.ancestry';
open OUT_ANCESTRY, "> $ancestry_file";
for my $pop (@pop_order){
my %repProp_ind;
    my $rep_K = $rep_pop_K{$pop};
    for my $ind (sort keys %{$pop_ind_line{$pop}}){
        my $infor = $qhash{$pop_ind_line{$pop}{$ind}};
        my $rep_K_infor = (split /\s+/, $infor)[$rep_K];
        $repProp_ind{$rep_K_infor}{$ind} = 1;
    }
    for my $repProp(reverse sort {$a<=>$b} keys %repProp_ind){
        for my $ind (sort keys %{$repProp_ind{$repProp}}){
            my $infor = $qhash{$pop_ind_line{$pop}{$ind}};
            print OUT_ANCESTRY "$ind\t$pop\t$infor\n";
        }
    }
}
close OUT_ANCESTRY;

#########################produce color outfile#############################
print "Producing .color file";
my $color_file = $outfile.'.color';
open COLOR_OUT, "> $color_file";
my @ancestry_number = split /\s+/, $qhash{'1'};
for my $i (0..$#ancestry_number){
    print COLOR_OUT "$colors[$i]\n";
}
close COLOR_OUT;

my $Rscript_file = $outfile.'.r';
my $pdf_file = $outfile.'.pdf';
my $png_file = $outfile.'.png';

print "Procuding Rscript.....\n";
open RSCRIPT, ">$Rscript_file";
print RSCRIPT "############################# define the sector function to draw sector ######################################\n";
print RSCRIPT "sector <- function (x0=0, y0=0, angle1, angle2, radius1, radius2, col, angleinc = 0.03)\n";
print RSCRIPT "{\nif (angle1 > angle2) {\ntemp <- angle1\nangle1 <- angle2\nangle2 <- temp\n}\nif (radius1 > radius2) {\ntemp <- radius1\nradius1 <- radius2\nradius2 <- temp\n}\n";
print RSCRIPT "###########use 4 points polygon to draw a sector########################\n";
print RSCRIPT "angles <- seq(angle1, angle2, by = angleinc)\nangles[length(angles)] <- angle2\nangles <- angles*(pi/180)\nxpos <- c(cos(angles) * radius1, cos(rev(angles)) * radius2) + x0\nypos <- c(sin(angles) * radius1, sin(rev(angles)) * radius2) + y0\npolygon(xpos, ypos, col = col, border = col)\n}\n";
print RSCRIPT "##############################################################################################################\n";
print RSCRIPT "args <- commandArgs(trailingOnly=T)\nfile <- args[1]\nout <- args[2]\nrmin <- 2\nrmax <- 1.85*rmin\namin <- -251\namax <- 91\nprgap <- 0.2\n";
print RSCRIPT "col <- c(";
for my $i (0..$#colors){
    if($i == 0){
        print RSCRIPT "'$colors[$i]'";
    }
    else{
        print RSCRIPT ",'$colors[$i]'";
    }
}
print RSCRIPT ")\nindcol <- 1\npopcol <- 2\nKstart <- 3\nKnum <- 1\n";

if(defined $opts{'f'}){
    print "defined -f png\n";
    print RSCRIPT "png(out, 1500, 1500)\npar('oma'=c(3,3,3,3))\npar('mar'=c(4, 4, 4, 4))\n";
}
else{
    print RSCRIPT "pdf(out, 45, 45)\npar('oma'=c(10,10,10,10))\npar('mar'=c(20, 20, 20, 20))\n";
}
print RSCRIPT "par('xpd'=TRUE)\nplot(0, 0, xlim=c(-rmax, rmax), ylim=c(-rmax, rmax), axes=F, ann=F, type='n')\nrstart <- rmin\nrlen <- (rmax-rmin)*(1-prgap)/Knum            ####### set the length unit\ndata <- read.table(file)\npopsize <- table(data[,popcol])\nnpop <- length(popsize)\nangelperpop <- (amax - amin)/npop\nangelperInd <- angelperpop/popsize\nangpre <- amin\n";

print RSCRIPT "#############################draw sectors##################################################################\n";
print RSCRIPT "for ( i in 1:nrow(data) ){                  ################ for each individual\n";
print RSCRIPT "ang1 <- angpre\nang2 <- ang1 + angelperInd[as.character(data[i,popcol])]\nrpre <- rstart\nfor ( j in Kstart:ncol(data) ){                  ################ for each K\nrpost <- rpre + rlen*data[i,j]\nsector(angle1=ang1, angle2=ang2, radius1=rpre, radius2=rpost, col=col[j-Kstart+1])   ###### draw sector for each K of each individual\nrpre <- rpost}\nangpre <- ang2\n}\nrstart <- rstart + (rmax-rmin)/Knum\n";

print RSCRIPT "lstart <- rmin\nlend <- rmax - prgap*(rmax-rmin)/Knum\n";
print RSCRIPT "#############get pop_number angles#################\nangles <- seq(amin, amax, length=npop+1)\n";
#####check if black_line exits############
if(defined $opts{'l'}){
    print "black lines separating populations are removed\n";
}
else{
print RSCRIPT "######################## draw lines for each pop #################################\n";
print RSCRIPT "for ( tmp in angles ){px <- c(lstart*cos(tmp*pi/180), lend*cos(tmp*pi/180))\npy <- c(lstart*sin(tmp*pi/180), lend*sin(tmp*pi/180))\nlines(px, py, col='black', cex=0.5)}\n";
}
print RSCRIPT "\n\n##########################write all population name#############################\n";
print RSCRIPT "unique <- as.character(unique(data[,popcol]))\nnpop <- length(unique)\nangelperpop <- (amax - amin)/npop\nangpre <- amin\n";
print RSCRIPT "for ( i in 1:npop ){\nang <- angpre + angelperpop/2\nxx <- lend*cos(ang*pi/180)\nyy <- lend*sin(ang*pi/180)\n";

if(defined $opts{'f'}){
    print RSCRIPT "text = unique[i];\ncex_no = 80/npop\nif(cex_no>1.6){cex_no=1.6}\nif(cex_no<1){cex_no=1}\n";
}
else{
    print RSCRIPT "text = unique[i];\ncex_no = 220/npop\nif(cex_no>6){cex_no=6}\nif(cex_no<1.5){cex_no=1.5}\n";
}
print RSCRIPT "if ( ang < -90 ){text(xx, yy, text, srt=180+ang, adj=c(1.1, 0.5), cex=cex_no, font=1, col=colors()[490])}\nelse {text(xx, yy, text, srt=ang, adj=c(-0.1, 0.5), cex=cex_no, font=1, col=colors()[490])\n}\n";


print RSCRIPT "angpre <- angpre + angelperpop\n}\n";


if(defined $opts{'t'}){
    print RSCRIPT "############################## draw piechart ######################################\n";
print RSCRIPT "orig <- data\n";
print RSCRIPT "rstart <- 0.0*rmin\nrend <- 0.4*rmin\nrlen <- rend-rstart\nx0 <- 0; y0 <- 0\namax <- 180; amin <- -180\ndata_pop <- orig[orig[,popcol]=='$target_pop',]\ndata_ind <- orig[orig[,indcol]=='$target_pop',]\ndata <- rbind(data_pop, data_ind)\napre <- amin\n";
print RSCRIPT "for ( j in Kstart:ncol(data) ){\nval <- mean(data[,j])\napost <- apre + (amax-amin)*val\nsector(x0=x0, y0=y0, angle1=apre, angle2=apost, radius1=rstart, radius2=rend, col=col[j-Kstart+1])\napre <- apost\nangles <- c(angles, apost)\n}\n";

if(defined $opts{'f'}){
    print RSCRIPT "text(0, 1, '$target_pop', cex=1, font=2, col=colors()[490])\n";
}
else{
    print RSCRIPT "text(0, 1, '$target_pop', cex=5, font=2, col=colors()[490])\n";
}
}
print RSCRIPT "dev.off()\n";
close RSCRIPT;

print "Drawing Figure.....\n";
if(defined $opts{'f'}){
    `Rscript $Rscript_file $ancestry_file $png_file`;
}
else{
    `Rscript $Rscript_file $ancestry_file $pdf_file`;
}

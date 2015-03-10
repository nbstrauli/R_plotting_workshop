setwd('~/Desktop/R_plotting_workshop')

plot.each.gene.expr = function(input_filepath){
	#This function plots each gene's expression trajectory over time as a semi-transparent line. 'input_filepath' is a tab-delimited file where each column represents a time-point, each row is a different gene, and each value in the matrix is that gene's expression level, at that time-point. The expression levels are normalized by time-point, so that each column sums to one.
	
	#First read in the data as a data frame
	t = read.table(input_filepath, header=TRUE)
	#remove last column of data because this is a summary statistic for each of the genes
	num_columns = length(t[1,])
	t = t[,-num_columns]
	
	#uncomment the line below if the data file has row identifiers as the 1st column (which is the case for the allele freq trajectory file)
	#t = t[,-1]
	
	#now find the maximum value in all the data to make appropriate limits for the plot's y-axis
	max_value = max(t)
	#also find number of time-points for x-axis limits
	num_tpoints = length(t[1,])
	
	#set filename and file-type (pdf in this case) for plot
	pdf(paste(input_filepath, '_line_plots.pdf', sep=''))
	
	#now use plot() function to plot first line of data, and set up the axes and such. The first two arguments give the x and y coordinants for each data-point, respectively.
	#see: https://stat.ethz.ch/R-manual/R-devel/library/graphics/html/plot.html
	#and https://stat.ethz.ch/R-manual/R-devel/library/graphics/html/par.html
	#for list of parameters that can be passed to plot()
	plot(x=0:(num_tpoints-1), y=t[1,], xlim=c(0,num_tpoints-1), ylim=c(0,max_value), main='Each V Gene Heavy Expression Trajectory', xlab='Time (days since vaccination)', ylab='Gene Expression Level', type='l', col=rgb(0,0,0,0.1))
	
	#now loop through all remaining rows (genes) in the data
	num_genes = length(t[,1])
	for (i in 2:num_genes){
		lines(x=0:(num_tpoints-1), y=t[i,], col=rgb(0,0,0,0.1))
	}
	
	#this will instruct R to close the file that the plot is being written
	dev.off()
}

plot.stacked.barplot = function(input_filepath){
	#This script uses the function barplot() to make a stacked bar plot of the data. Here, the overall length of a bar for a time-point corresponds to the cumulative gene expression at that time. Each color within a bar corresponds to a unique gene, so that the length of an individual colored segment corresponds to that gene's expression level
	#read in data
	t = read.table(input_filepath, header=TRUE)
	#convert to matrix because 'barplot()' only takes vectors of matrices as input.
	m = as.matrix(t, mode='numeric')
	
	#make different colors for different genes using rainbow() function
	#see: https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/palettes.html for more info
	#note that this only generates 10 colors. So the color pattern repeats every 10 genes, but this should be sufficient to distinguish between genes in the plot.
	colors = rainbow(10)
	
	#remove last column because it is not a time-point
	num_columns = length(m[1,])
	m = m[,-num_columns]
	num_tpoints = length(m[1,])
	
	#now make barplot. See: https://stat.ethz.ch/R-manual/R-devel/library/graphics/html/barplot.html for more info
	pdf(paste(input_filepath, '_stacked_barplot.pdf', sep=''))
	barplot(m, col=colors, names.arg=0:(num_tpoints-1), main='Each V Gene Heavy Expression Level', xlab='Time (days since vaccination)', ylab='Cumulative Expression Level')
	dev.off()
}

parse.data.for.ggplot = function(input_filepath){
	#This script will parse the data with in input filepath (which is essentially a matrix) into a large data frame (which is the proper input for ggplot)
	#read in data and then turn it into a matrix
	t = read.table(input_filepath, header=TRUE)
	num_columns = length(t[1,])
	m = as.matrix(t[,-num_columns], mode='numeric')
	
	#get dimensions of the matrix
	num_genes = length(m[,1])
	num_tpoints = length(m[1,])
	
	#now we need to cycle through the elements of the matrix and append to our new data frame each time
	#each of the empty variables below will become the columns to our data frame
	expr_level = c()
	tpoints = c()
	gene_ids = c()
	
	#1st cycle through each row (i.e. gene of the matrix)
	for (i in 1:num_genes){
		#assign a random ID for this gene (in this case a random number b/t 0 and 100)
		gene_id = runif(1, 0, 100)
		gene_id = as.character(gene_id)
		tpoint = 0
		
		#now cycle through each element of the given row (i.e. time-point)
		for (j in m[i,]){
			#append each of our variables (i.e. soon-to-be columns in our data frame)
			expr_level = append(expr_level, j)
			tpoints = append(tpoints, tpoint)
			gene_ids = append(gene_ids, gene_id)
			
			#update the timepoint each with each cycle
			tpoint = tpoint + 1
		}
	}
	
	#now combine each of the variables so that they are columns in a data frame
	data = data.frame(expr_level, tpoints, gene_ids)
	
	#write the data frame to file using write.table()
	#see: https://stat.ethz.ch/R-manual/R-devel/library/utils/html/write.table.html for more info
	data_filepath = paste(input_filepath, '_ggplot_friendly_dataframe', sep='')
	write.table(data, data_filepath, row.names=FALSE, quote=FALSE)
	
	#return the data frame as well
	return (data)
}

plot_stacked_area = function(input_filepath){
	#This script makes a stacked area chart using ggplot2's geom_area() function. These plots essentially convey the same information as a stacked bar plot but are way prettier, and thus more interpretable. Here, each color corresponds to a gene, and the width of that color at a given time-point gives it expression level for that time.

	#1st the data needs to be parsed in a way that is amenable for ggplot plotting
	data_frame = parse.data.for.ggplot(input_filepath)

	#must load ggplot may need to use install.packages('ggplot2') if ggplot has not been installed yet
	library(ggplot2)

	#This tells ggplot how to interperate the data, i.e. which columns in the data frame give the x coordinants, which give the y, etc...
	p = ggplot(data_frame, aes(x=tpoints, y=expr_level, fill=gene_ids))
	#this tells ggplot what type of plot you would like to make
	p = p + geom_area(position='stack')
	#this instructs ggplot to not include a legend
	p = p + guides(fill=FALSE)
	#this gives the title and axis labels for the plot
	p = p + ggtitle('Each V Gene Expression Trajectories') + xlab('Time (days since vaccination)') + ylab('Cumulative Expression Level')

	#now save the plot
	pdf(paste(input_filepath, '_stacked_area_plot.pdf', sep=''))
	plot(p)
	dev.off()
}

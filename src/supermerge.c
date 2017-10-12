/* Copyright (C) 2016-2017 Bohdan Khomtchouk and Derek Van Booven */

/* This file is part of SUPERmerge. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int debug = 0;

main(int argc, char **argv)
{
	extern char *optarg;
	extern int optind;
	int c, err = 0;
	char *length = "default_depth";
	char *length2 = "default_interval_distance";
	static char usage[] = "usage: %s [-d] depth_cutoff [-i] interval_distance \n";
	char *gtf_file;
	char system_command[1000];
	int return_value;

	while ((c = getopt(argc, argv, "d:i:g:")) != -1)
		switch (c) {
			case 'd':
				length = optarg;
				break;
			case 'i':
				length2 = optarg;
				break;
			case 'g':
				gtf_file = optarg;
				break;
		}
	if (strcmp(length,"default_distance") == 0) {
		printf("no parameter was specified so default to 20\n");
	}
	else {
		printf("depth input was %s\n", length);
	}

	/* Input was the bam file so now take the bam file and pass it to the coverage function for chopping */
	get_coverage(argv[optind], atoi(length));
	get_intervals(argv[optind], atoi(length2), atoi(length));
	get_annotation(argv[optind], gtf_file);
	

	// Cleanup the temporary files.
	sprintf(system_command, "cut -f5 %s.tmp3 > percentages.stats", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "grep protein_coding %s.tmp3 | wc -l >> annotation.stats", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "grep lincRNA %s.tmp3 | wc -l >> annotation.stats", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "grep miR %s.tmp3 | wc -l >> annotation.stats", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "grep antisense %s.tmp3 | wc -l >> annotation.stats", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "grep pseudogene %s.tmp3 | wc -l >> annotation.stats", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "cut -f6 %s.tmp2 > length.stats", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "cut -f7 %s.tmp2 > libsize.stats", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "rm %s.tmp2", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "rm %s.tmp", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "rm tmp.gtf", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "mv %s.tmp3 %s.results", argv[optind], argv[optind]);
	return_value = system(system_command); 

	sprintf(system_command, "echo \"x<-scan(\\\"percentages.stats\\\")\" > graphs.R", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "echo \"hist(x, main=NULL, xlab=\\\"Percentage\\\")\" >> graphs.R", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "echo \"y<-scan(\\\"annotation.stats\\\")\" >> graphs.R", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "echo \"z<-c(\\\"protein_coding\\\", \\\"lincRNA\\\", \\\"miR\\\", \\\"antisense\\\", \\\"pseudogene\\\")\" >> graphs.R", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "echo \"barplot(y, names.arg=z, main=NULL, ylab=\\\"Frequency\\\")\" >> graphs.R", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "echo \"w<-scan(\\\"length.stats\\\")\" >> graphs.R", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "echo \"hist(log(w), main=NULL, xlab=\\\"log(coverage island length)\\\")\" >> graphs.R", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "echo \"t<-scan(\\\"libsize.stats\\\")\" >> graphs.R", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "echo \"boxplot(log(t), main=NULL, ylab=\\\"log(library size)\\\")\" >> graphs.R", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "echo \"dev.off()\" >> graphs.R", argv[optind]);
	return_value = system(system_command);

	sprintf(system_command, "R CMD BATCH graphs.R", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "mv Rplots.pdf %s.pdf", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "rm graphs.R", argv[optind]);
	return_value = system(system_command); 
	sprintf(system_command, "rm *stats", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "rm *.Rout", argv[optind]);
	return_value = system(system_command);
	sprintf(system_command, "rm *.pileup", argv[optind]); 
	return_value = system(system_command);

}


int get_coverage(char *input_filename, int depth)
{

/* This function will have a system command to run samtools and reduce mpileup file into 3 columns for merging.*/
/* If necessary we can put the samtools ripped code here to produce our own mpileup function for depth.*/

//	printf("system command is samtools mpileup %s | cut -f1,2,4 | awk '{ if ($3 > %d) print }'", input_filename, depth);
	int return_value;
	char system_command[1000];
	//sprintf(system_command, "samtools mpileup %s | cut -f1,2,4 | awk '{ if ($3 > %d) print }' > %s.tmp", input_filename, depth, input_filename);
	sprintf(system_command, "samtools mpileup %s > %s.pileup", input_filename, input_filename);
	return_value = system(system_command);
	sprintf(system_command, "cut -f1,2,4 %s.pileup | awk '{ if ($3 > %d) print }' > %s.tmp", input_filename, depth, input_filename);
	return_value = system(system_command);

	return return_value;

}

int get_intervals(char *input_filename, int intervals, int depth3)
{
/*  This function takes the output from the coverage function and creates the specified intervals. */

	//open file and now expand the intervals by desiered distance
	
	FILE *depthfile, *tmp2file, *finalfile;

	char system_file_tmp[1000];
	char file_line[100000], chromosome[10], position[32], depth2[5];
	char chromosome_current[10];
	char *token;
	long start, stop, positions, running_total;

	sprintf(system_file_tmp, "%s.tmp", input_filename);
	depthfile = fopen(system_file_tmp, "r");
	sprintf(chromosome_current, "0x0");
	start = 0;
	stop = 0;
	positions = 1;
	running_total = 0;
	if (depthfile == NULL) {
		printf("problem opening temp depth file\n");
	} else {
		sprintf(system_file_tmp, "%s.tmp2", input_filename);
		tmp2file = fopen(system_file_tmp, "w");
		if (fgets(file_line, 100000, depthfile) == NULL) {
			printf("we have a null situation, no reported intervals\n");
			return 5;
		}
		sscanf(file_line,"%s%s%s", chromosome, position, depth2);
		sprintf(chromosome_current, "%s", chromosome);
		start = atol(position) - intervals;
		stop = atol(position) + intervals;
		running_total = atol(depth2);
		while (fgets(file_line, 100000, depthfile) != NULL) {
			sscanf(file_line, "%s%s%s", chromosome, position, depth2);
			// check if chromosome is the same as the previous line chromosome, or if we have a new chromosome.
			if (strcmp(chromosome, chromosome_current) != 0) {
				// new chromosome so report what we see in the interval and reset the values.
				fprintf(tmp2file, "%s\t%d\t%d\t%d\t%f\t%d\n", chromosome_current, start, stop, positions, (float)positions / (float)(1 + stop - start), running_total);
				sprintf(chromosome_current, "%s", chromosome);
				start = atol(position) - intervals;
				stop = atol(position) + intervals;
				positions = 1;
				running_total = atol(depth2);
			}
			else {
				//we are in the same chromosome so check the start position for the inclusion.
				if (((atol(position)-intervals) >= start) && ((atol(position) - intervals) <= stop)) {
					//we have inclusion so set new stop to the current stop
					stop = atol(position) + intervals;
					positions = positions + 1;
					// Increment running_total appropriately.
					running_total = running_total + atol(depth2);
				}
				else {
					//we do not have inclusion so print the current interval and reset starting values to line of file
					fprintf(tmp2file, "%s\t%d\t%d\t%d\t%f\t%d\t%d\n", chromosome_current, start, stop, positions, (float)positions / (float)(1 + stop - start), (1 + stop - start), running_total);
					sprintf(chromosome_current, "%s", chromosome);
					start = atol(position) - intervals;
					stop = atol(position) + intervals;
					positions = 1;
					running_total = atol(depth2);
				}
			}
		}
		// don't forget to print the last interval
		fprintf(tmp2file, "%s\t%d\t%d\t%d\t%f\t%d\t%d\n", chromosome_current, start, stop, positions, (float)positions / (float)(1 + stop - start), ( 1 + stop - start), running_total);
		positions = 1;
		running_total = 0;
	}


	
	//cleanup
	fclose(depthfile);
	fclose(tmp2file);
	return 5;
}

int get_annotation(char *input_filename, char *gtf)
{

	char system_command[1000];

	sprintf(system_command,"awk -F \"\\t\" '$3 == \"gene\" { print $1 \"\\t\" $4 \"\\t\" $5 \"\\t\" $9 }' %s > tmp.gtf", gtf);

	system(system_command);

	// remove the chr if they are in the file
	sprintf(system_command, "sed -i 's/^chr//g' tmp.gtf");
	system(system_command);

	printf("we have a gtf file with the name of %s\n", gtf);

	char varstr[50000];
	char vcfstr[50000];
	char varstr1[50000];
	char vcfstr1[50000];
	char tempstring[50000];
	char col1[100];
	char col2[100];
	char col3[100];
	char col4[100];
	char col5[100];

	char * pvarcol1;
	char * pvarcol2;
	char * pvarcol3;
	char * pvarcol4;
	char * pvarcol5;
	char * pvarcol6;
	char * pvcfcol1_1;
	char * pvcfcol2_2;
	char * pvcfcol3_3;
	char * pvcfcol4_4;
	char * pvcfcol5_5;
	char * pvcfcol1;
	char * pvcfcol2;
	char * pvcfcol3;
	char * pvcfcol4;
	char * pvcfcol5;
	char * end_str;
	char * end_str2;
	char * commentptr;
	FILE * vcffp;
	FILE * varfp;
	FILE * outfp;
	long i;
	long j;	
	long j_1;
	long k;
	long l;
	long l_1;
	long m;
	long n;
	long n_1;
	long v1;
	long v2;
	long v3;
	long v4;
	int zero_count;
	int pflag;
	int wflag;

	wflag = 1;
	zero_count = 0;
	pflag = 0;
	/* open up all the files and make sure they are open, else return value 1 */
	sprintf(system_command, "sed -i 's/^chr//g' %s.tmp2", input_filename);
	system(system_command);
	sprintf(system_command, "sed -i 's/^X/23/g' %s.tmp2", input_filename);
	system(system_command);
	sprintf(system_command, "sed -i 's/^Y/24/g' %s.tmp2", input_filename);
	system(system_command);
	sprintf(system_command, "sed -i 's/^M/25/g' %s.tmp2", input_filename);
	system(system_command);
	sprintf(system_command, "sed -i 's/^MT/25/g' %s.tmp2", input_filename);
	system(system_command);
	sprintf(system_command, "sort -k1,1n -k2,2n %s.tmp2 > tmp.tmp2", input_filename);
	system(system_command);
	sprintf(system_command, "mv tmp.tmp2 %s.tmp2", input_filename);
	system(system_command);
	sprintf(tempstring,"%s.tmp2", input_filename);

        varfp = fopen(tempstring, "r");
	if (varfp == NULL)
	{
		printf("Problem opening variant file\n");
		return 1;
	}
        vcffp = fopen("tmp.gtf", "r");
	if (vcffp == NULL)
	{
		printf("Problem opening vcf file\n");
		return 1;
	}

	sprintf(tempstring,"%s.tmp3", input_filename);
	outfp = fopen(tempstring,"w");
	if (outfp == NULL)
	{
		printf("Problem opening output file\n");
	}
	// Get the first lines of both files and see what's up.

	fgets(varstr, 50000, varfp);
	strcpy(varstr1,varstr);
        pvarcol1 = strtok_r(varstr,"\t\n",&end_str);
	k = strtol(pvarcol1,NULL,0);
        pvarcol2 = strtok_r(NULL,"\t\n",&end_str);
        i = strtol(pvarcol2,NULL,0);
        pvarcol3 = strtok_r(NULL,"\t\n",&end_str);
        m = strtol(pvarcol3,NULL,0);
        pvarcol4 = strtok_r(NULL,"\t\n",&end_str);
        pvarcol5 = strtok_r(NULL,"\t\n",&end_str);
        pvarcol6 = strtok_r(NULL,"\t\n",&end_str);


	fgets(vcfstr, 10000, vcffp);
	strcpy(vcfstr1,vcfstr);
        pvcfcol1 = strtok_r(vcfstr,"\t",&end_str2);
	n = strtol(pvarcol1,NULL,0);
        pvcfcol2 = strtok_r(NULL,"\t",&end_str2);
	j = strtol(pvcfcol2,NULL,0);
        pvcfcol3 = strtok_r(NULL,"\t",&end_str2);
        l = strtol(pvcfcol3,NULL,0);
        pvcfcol4 = strtok_r(NULL,"\t",&end_str2);
        pvcfcol5 = strtok_r(NULL,"\t\n",&end_str2);

	pvcfcol1_1 = pvcfcol1;
	pvcfcol2_2 = pvcfcol2;
	pvcfcol3_3 = pvcfcol3;
	pvcfcol4_4 = pvcfcol4;
	pvcfcol5_5 = pvcfcol5;	
	strcpy(col1,pvcfcol4);
	n_1 = n;
	j_1 = j;
	l_1 = l;

	v1 = i - j;
        v2 = i - l;
        v3 = m - j;
        v4 = m - l;



	while(wflag > 0)
	{

		// check chromosomes.  If gtf is on next chromosome then move var file appropriately.
		
		if (k < n) {
			fprintf(outfp,"%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%d\n", pvarcol1, pvarcol2, pvarcol3, pvarcol4, pvarcol5, pvarcol6, n_1, j_1, l_1, col1, abs((i - l_1))+1);
			pflag = 1;
		}
		else if (k > n) {
			// Increase the gtf file since there are remaining chromosomes.
			pflag = 2;

		}
		else {
			if ((v1 > 0) && (v2 > 0) && (v3 > 0) && (v4 > 0)) {
				// if the peak in the gtf is all positive then we need to update the previous pointer to the current gtf peak

				pvcfcol1_1 = pvcfcol1;
				pvcfcol2_2 = pvcfcol2;
				pvcfcol3_3 = pvcfcol3;
				pvcfcol4_4 = pvcfcol4;
				
				strcpy(col1,pvcfcol4);
				n_1 = n;
				j_1 = j;
				l_1 = l;
				// set pflag to increase the gtf file
				pflag = 2;
			}
			else if ((v1 < 0) && (v2 < 0) && (v3 < 0) && (v4 < 0)) {
				// if all is negative then we need to which which is the closest peak and then move the peak file to the next flag.
				if ((j - m) < (i - l_1)) {
					// the previous pointer is closre to the distance so we need to put the distance there
					fprintf(outfp,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\n", pvarcol1, pvarcol2, pvarcol3, pvarcol4, pvarcol5, pvarcol6, pvcfcol1, pvcfcol2, pvcfcol3, pvcfcol4, abs((j - m))+1);
				} else {
					if (n_1 == n) {
						fprintf(outfp,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%d\n", pvarcol1, pvarcol2, pvarcol3, pvarcol4, pvarcol5, pvarcol6, pvcfcol1, j_1, l_1, col1, abs((i - l_1))+1);
					}
					else {
						fprintf(outfp,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\n", pvarcol1, pvarcol2, pvarcol3, pvarcol4, pvarcol5, pvarcol6, pvcfcol1, pvcfcol2, pvcfcol3, pvcfcol4, abs((j - m))+1);
					}

					
				}
				// set the pflag to represent a movement in the peak file
				pflag = 1;
			}
			else {
			/* add to the running total of the overlap and move ot the next query */
				zero_count = zero_count + 1;
				fprintf(outfp,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t0\n", pvarcol1, pvarcol2, pvarcol3, pvarcol4, pvarcol5, pvarcol6, pvcfcol1, pvcfcol2, pvcfcol3, pvcfcol4);
				pflag = 1;
			}
		}

		switch(pflag) {
			case 1 :
				// increase the gtf file
				if (fgets(varstr, 50000, varfp) != NULL) {
                                        strcpy(varstr1,varstr);
                                        pvarcol1 = strtok_r(varstr,"\t\n",&end_str);
                                        k = strtol(pvarcol1,NULL,0);
                                        pvarcol2 = strtok_r(NULL,"\t\n",&end_str);
                                        i = strtol(pvarcol2,NULL,0);
                                        pvarcol3 = strtok_r(NULL,"\t\n",&end_str);
                                        m = strtol(pvarcol3,NULL,0);
                                        pvarcol4 = strtok_r(NULL,"\t\n",&end_str);
                                        pvarcol5 = strtok_r(NULL,"\t\n",&end_str);
                                        pvarcol6 = strtok_r(NULL,"\t\n",&end_str);
                                } else { wflag = 0; }; 
				break;
			case 2 :
				// increase the peak file
				if (fgets(vcfstr, 10000, vcffp) != NULL) {
                                        strcpy(vcfstr1,vcfstr);
                                        pvcfcol1 = strtok_r(vcfstr,"\t",&end_str2);
                                        n = strtol(pvcfcol1,NULL,0);
                                        pvcfcol2 = strtok_r(NULL,"\t",&end_str2);
                                        j = strtol(pvcfcol2,NULL,0);
                                        pvcfcol3 = strtok_r(NULL,"\t",&end_str2);
                                        l = strtol(pvcfcol3,NULL,0);
                                        pvcfcol4 = strtok_r(NULL,"\t\n",&end_str2);
                                        pvcfcol5 = strtok_r(NULL,"\t\n",&end_str2);
                                } else { wflag = 0; }
				break;
		}

		v1 = i - j;
		v2 = i - l;
		v3 = m - j;
		v4 = m - l;
	}
	fclose(varfp);
	fclose(vcffp);
	return 0;
}


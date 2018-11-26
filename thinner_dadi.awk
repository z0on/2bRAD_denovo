# thins dadi formatted file to specified interval between markers 
# dadi file must be sorted by chromosome (contig) and position
# default interval: 250
# usage:
# awk -f thinner_dadi.awk -v interval=500 my_dadi.data >dadi_thinned.data
BEGIN {
	CHR=""
	POS=0
#	print("interval:",interval)
	if (length(interval)==0) 
	{ 
		interval=250
		print("using default interval:",interval)
	}
}

{
# print("----------\nold:",CHR,POS,"\nnew:",$(NF-1),$(NF))
if ($(NF-1)!=CHR) 
	{
	CHR=$(NF-1)
	POS=$(NF)
	print
	}
else 
	{
	srand($(NF))
	roll=rand()
	test=(($(NF)-POS)/interval)^3
#	print("roll:",roll,"\ntest:",test)
	if(roll<test)
		{	
		POS=$(NF)
		print 
		}
#	else 
#		{
#			print ("\tskipping ", $(NF-1),$(NF))
#		}
	}
}

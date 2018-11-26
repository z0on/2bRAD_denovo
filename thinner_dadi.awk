# thins dadi formatted file to specified interval between markers *on average*
# different runs produce slightly different thinning 
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
	srand()
}
{
if ($(NF-1)!=CHR) 
	{
	CHR=$(NF-1)
	POS=$(NF)
	print
	}
else 
	{
	if(rand()<(($(NF)-POS)/interval)^3)
		{	
		POS=$(NF)
		print 
		}
#	else 
#		{
#			print ("skipping ", $(NF-1),$(NF))
#		}
	}
}

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
	if(rand()<($(NF)-POS)/(2*interval))
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

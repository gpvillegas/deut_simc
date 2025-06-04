! HISTOGRAMS.INC
! Storage space for the SIMULATE histogram arrays

! Record structures that are only needed for histograms

	integer		nHbins
	parameter	(nHbins=50)

	type hist_entry
    	   sequence
	  real*8	bin,min
	  real*8	buf(nHbins)
	end type

	type hist_arm
	    sequence
	    type(hist_entry):: delta, yptar, xptar
        end type

	type hist_arm2
            sequence
	    type(hist_entry):: delta, yptar, xptar
        end type

	type hist_double_arm 
          sequence
	  type(hist_arm):: e
	  type(hist_arm2)::p
	  type(hist_entry):: Em, Pm
	end type

	type histograms
          sequence	
	  type(hist_double_arm):: RECON, gen, geni
	end type

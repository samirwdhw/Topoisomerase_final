c	n enz, adjustive force, n defined

	parameter (Vm=3.6, Km=270, T=300, Vo=3.38, ATP=2000, d=1, nn=1, nEnzyme=500)
	parameter (maxconf= 4000, dt=0.001, bp=30000, maxcatenation=11*(bp/100))
	parameter (ap=0.241432)				  !ap is angie's parameter for KbT and pico newton settlement
	integer chro(bp), timep(maxconf), timesum
	double precision it,turnover,DNAextension
	open(unit=20, file='tracingchroDP_OnCen_Ref2_500e30kbp03.dat')
	open(unit=10, file='chroDP_OnCen_Ref2_500e30kbp2000a03.dat')	
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		
c		avo=6.0221415E23
c		pico=6.0221415E11
c		n=nn*pico
		
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
		
		do 1 iforce= 0, 3 ! 10			!force loop
		fef=iforce		
		timesum=0
		count=0
c		print*, "Force=", iforce
		
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		do 2 iconf=1,maxconf		!config loop
		it=0			
				
				do 11 icro=1,bp	!initialising to zero
				chro(icro)=0
11				continue
				
			        icat2=maxcatenation

      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc			        
			
			        icat2=maxcatenation
			
				do 22 iiins= 1,maxcatenation 
					
33						iins = (4*bp*0.1) + rand()* bp*2*0.1
							if (chro(iins).eq.0) then
							chro(iins)=1
					
							else if (chro(iins).eq.1) then
							goto 33
							endif					!insertion of catenations

						
22				continue
				


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
1000	it=it+1
	
		do 37 ien= 1, nEnzyme
	
				iengo=1+bp*rand()
					
				if ( chro(iengo).eq.1 ) then
				z=(abs((igo-(bp/2))/bp))
				feff= z*(fef)/icat2
					
				Vf= Vo*exp(-feff*ap)*(ATP)/(ATP+Km)
				prob=Vf*dt
					
				r1=rand()
				
				if (r1.le.prob) then
				chro(iengo)=0
				count=count+1
				icat2=icat2-1
c				print*, "t=",it*dt, " icat2=",icat2
				endif
				endif
						
37		continue

	if(icat2.ne.0)then
		goto 1000					!time loop ends
	else
		goto 1001
	end if


1001			timep(iconf)= it		!maximum time

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c			print*, icat2
			if(mod(iconf,2).eq.0) print*,"F=",iforce," C=",iconf," t=",it,
     1                  " icat2=",icat2," count=",count			
2			continue				!config loop ends

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


1111		do 19 iconf=1,maxconf
		timesum=timesum+(timep(iconf)*dt)
19		continue
		avgtime=timesum/maxconf

		write(20,*) 'Force loop= ',iforce
		
		turnover= count/(timesum*nEnzyme)
		DNAextension=turnover*90*maxcatenation				! 90 nm is the extension after resolving 1 catenation
		DNAextperunittime=DNAextension/avgtime				! DNA Extension per unit time 

		print*, fef, turnover, DNAextperunittime, (timesum)/maxconf        !icount, timesum
		write (10,*) fef, turnover, DNAextperunittime, (timesum)/maxconf     !,icount, timesum
		
1		continue 				!force loop ends


	stop
	end			
	

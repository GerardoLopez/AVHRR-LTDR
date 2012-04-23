#!/bin/csh -f

set day = $1
set month = $2
set year = $3

#--- Function to transform DDD YYYY to to DD.MM.YYYY ej. 2008.01.01 ---#
if ( $year == 1980 || $year == 1984 || $year == 1988 || $year == 1992 || $year == 1996 || $year == 2000 || $year == 2004 || $year == 2008 || $year == 2012 ) then
	set months_days = (31 29 31 30 31 30 31 31 30 31 30 31)
else
	set months_days = (31 28 31 30 31 30 31 31 30 31 30 31)
endif

set days_counter = 0
if ($month == "01") then
	set julian_day = `echo $day + 0 | bc -l`
else
	set month = `echo $month - 1 | bc -l`
	foreach i (`seq 1 1 $month`)
		@ days_counter = $days_counter + $months_days[$i]
	end	
	set julian_day = `echo $days_counter + $day | bc -l`
endif

echo $julian_day

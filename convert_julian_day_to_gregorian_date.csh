#!/bin/csh -f

set julian_day = $1
set year = $2

#--- Function to transform DDD YYYY to to DD.MM.YYYY ej. 2008.01.01 ---#
if ( $year == 1980 || $year == 1984 || $year == 1988 || $year == 1992 || $year == 1996 || $year == 2000 || $year == 2004 || $year == 2008 || $year == 2012 ) then
	set months_days = (31 29 31 30 31 30 31 31 30 31 30 31)
else
	set months_days = (31 28 31 30 31 30 31 31 30 31 30 31)
endif

set counter = 1
set flag = 0

set days_counter = $months_days[1]

while ( $flag == 0)
	if ($julian_day > $days_counter) then
		@ counter ++
		@ days_counter = $days_counter + $months_days[$counter]
	else
		@ day = $months_days[$counter] - $days_counter + $julian_day
		@ month = $counter
		set flag = 1
	endif
end

set day = `echo $day | awk '{if (length($1)==1) $1="0"$1;} {print $1}'`
set month = `echo $month | awk '{if (length($1)==1) $1="0"$1;} {print $1}'`

echo $year$month$day
